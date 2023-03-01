import os

import pandas as pd
import numpy as np
from cobra.io import read_sbml_model
from cobra.sampling import ACHRSampler
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
from scipy import stats
from scipy.stats import hypergeom, ks_2samp
from statsmodels.stats import multitest
import sys
# sys.path.insert(0, r"C:\Users\Bisbii\PythonProjects\ExpAlgae\src")
sys.path.insert(0, "/home/src/")
from ExpAlgae.model.COBRAmodel import MyModel


def achr_sample(filename, biomass_reaction):
    print('Sampling model: ', filename)
    folder = '/'.join(filename.split("/")[:-1])
    model_to_sample = read_sbml_model(filename)
    model_to_sample.objective = biomass_reaction
    # initial_solution = model_to_sample.optimize().objective_value
    model_to_sample.reactions.get_by_id(biomass_reaction).lower_bound = 0.15
    sampler = ACHRSampler(model_to_sample, thinning=100, seed=42)
    samples = sampler.sample(10000)
    result_filename = f"{folder}/{filename.split('/')[-1].split('.xml')[0]}_ACHR_samples.csv"
    samples.to_csv(result_filename, index=False)


def load_results(file_names):
    data_list = []
    for filename in file_names:
        dataframe = pd.read_csv(filename)
        if "Unnamed: 0" in dataframe.columns:
            dataframe.drop(columns=["Unnamed: 0"], inplace=True)
        # dataframe = dataframe * dataframe[biomass_map[filename.split(".")[0]]]
        data_list.append(dataframe)
    reactions = list(set([reac for dataframe in data_list for reac in dataframe.columns]))
    for reaction in reactions:
        for index, value in enumerate(data_list):
            if reaction not in value.columns:
                data_list[index][reaction] = 0
    results = {}
    for reaction in reactions:
        data = [value for dataframe in data_list for value in [dataframe[reaction]]]
        if not np.allclose(data[0], data[1]):
                # and not np.all_close(data[1], data[2]) and  not np.all_close(data[2], data[3]):
            stat, p_value = stats.kstest(data[0], data[1])
            print(f"Reaction: {reaction}")
            print("Statistic: ", stat)
            print("p-value: ", p_value)
            results[reaction] = [stat, p_value]
    results_df = pd.DataFrame.from_dict(results, orient="index", columns=["Statistic", "p-value"])
    results_df.to_csv("ACHR_results.csv")

def kstest(samples_healthy: pd.DataFrame, samples_infected: pd.DataFrame, dataset_name: str):
    """
    Calculate the K-S test to detect significantly altered reactions fluxes.
    Results are saved in a csv file.

    Parameters
    ----------
    samples_healthy: pd.DataFrame
        The samples of the healthy tissue.
    samples_infected: pd.DataFrame
        The samples of the infected tissue.
    dataset_name: str
        The name of the dataset.
    Returns
    -------
    pd.DataFrame: The results of the K-S test for each reaction.

    """
    rxns1 = set(samples_infected.columns)
    rxns2 = set(samples_healthy.columns)

    rxn_c = rxns1.intersection(rxns2)

    pvals = []
    rxnid = []
    fc = []

    for rxn in rxn_c:
        data1 = samples_infected[rxn].round(decimals=4)
        data2 = samples_healthy[rxn].round(decimals=4)

        data1 = data1.sample(n=1000)
        data2 = data2.sample(n=1000)

        if (data1.std() != 0 and data1.mean() != 0) or (data2.std() != 0 and data2.mean() != 0):
            kstat, pval = ks_2samp(data1, data2)

            foldc = (data1.mean() - data2.mean()) / abs(data1.mean() + data2.mean())

            pvals.append(pval)
            rxnid.append(rxn)
            fc.append(foldc)

    data_mwu = pd.DataFrame({'Reaction': rxnid, 'Pvalue': pvals})
    data_mwu = data_mwu.set_index('Reaction')

    reject, padj, _, _ = multitest.multipletests(data_mwu['Pvalue'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    data_mwu['Padj'] = padj
    data_mwu['Reject'] = reject
    data_mwu['FC'] = fc

    data_sig_fc = data_mwu.loc[(abs(data_mwu['FC']) > 0.82) & (data_mwu['Padj'] < 0.05), :]

    rxns1 = set(samples_infected.columns)
    rxns2 = set(samples_healthy.columns)

    rxn_in1 = rxns1.difference(rxns2)
    rxn_in2 = rxns2.difference(rxns1)

    sigs = Parallel(n_jobs=8)(delayed(bootstrap_ci)(samples_infected[rx]) for rx in rxn_in1)
    act = [sigs[i][0] for i in range(len(sigs)) if sigs[i][1] == 1]
    sigs2 = Parallel(n_jobs=8)(delayed(bootstrap_ci)(samples_healthy[rx]) for rx in rxn_in2)
    rep = [sigs2[i][0] for i in range(len(sigs2)) if sigs2[i][1] == 1]

    df_abs = pd.DataFrame({'Reaction': act + rep, 'Padj_bootstrap': np.zeros(len(act + rep))})
    df_abs = df_abs.set_index('Reaction')
    # data_return = data_sig_fc + df_abs
    data_return = pd.concat([data_sig_fc, df_abs])
    file = f"{dataset_name}_results.csv"
    data_return.to_csv(file)

    return data_return.index.to_list()


def bootstrap_ci(rxn):
    """
    Calculate the confidence interval of a reaction

    Parameters
    ----------
    rxn

    Returns
    -------
    int: 1 if the reaction is significantly different from the mean, 0 otherwise.

    """
    bsci = []

    for i in range(1000):
        bt_samp = rxn.sample(1000, replace=True)
        bsci.append(bt_samp.mean())

    ci_low = np.percentile(bsci, 2.5)
    ci_high = np.percentile(bsci, 97.5)

    if ci_low > 0 or ci_high < 0:
        return rxn.name, 1
    else:
        return rxn.name, 0


def pathway_enrichment(rxnlist: list, dataset_name: str):
    """
    Maps significantly altered reactions to pathways using the subsystems from the HumanGEM model.
    Results are saved in csv and jpg files.

    Parameters
    ----------
    rxnlist: list
        The list of reactions to be mapped to pathways.
    dataset_name: str
        The name of the dataset.
    """

    dataset = pd.read_csv("pathways_map.csv")

    listrxn_size = []
    set_size = []

    d = [g for g in rxnlist]

    for col in dataset.columns:
        dataframe = pd.DataFrame({'Reaction': dataset[col]})

        out = []

        for reac in dataframe['Reaction']:
            if reac in rxnlist:
                out.append(reac)
                if reac in d:
                    d.remove(reac)

        listrxn_size.append(len(out))
        set_size.append(len(dataset[col].dropna()))

    hyperdata = pd.DataFrame({'Pathways': dataset.columns, 'ListReactions': listrxn_size, 'SetSize': set_size})

    hits = hyperdata['ListReactions']
    pool = hyperdata['SetSize']

    allrxns = hyperdata['SetSize'].sum()
    targetrxns = hyperdata['ListReactions'].sum()

    pval_list = []

    for h, p in zip(hits, pool):
        rv = hypergeom(allrxns - p, p, targetrxns)

        pval = rv.pmf(h)

        pval_list.append(pval)

    hyperdata['P-value'] = pval_list

    reject, padj, _, _ = multitest.multipletests(hyperdata['P-value'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

    hyperdata['P-value_adj'] = padj
    hyperdata['Reject'] = reject

    hyperdata_sig = hyperdata[(hyperdata['Reject']) & (hyperdata['ListReactions'] != 0)]

    hyperdata_sorted = hyperdata_sig.sort_values(by='P-value_adj', ascending=False)
    to_ignore = ['Transporters pathway', 'Drains pathway', 'Biosynthesis of cofactors', 'Microbial metabolism in diverse environments']
    hyperdata_sorted = hyperdata_sorted[~hyperdata_sorted['Pathways'].isin(to_ignore)]

    plt.figure(figsize=(12, 10))

    sc = plt.scatter(hyperdata_sorted['P-value_adj'], np.arange(0, len(hyperdata_sorted['Pathways'])), s=hyperdata_sorted['ListReactions'], color=(0.9, 0.3, 0.1, 0.9))

    plt.xlabel('Adjusted p-value')

    plt.yticks(np.arange(0, len(hyperdata_sorted['Pathways'])), labels=hyperdata_sorted['Pathways'])

    handles, labels = sc.legend_elements(prop="sizes", alpha=0.8)

    plt.legend(handles, labels, bbox_to_anchor=(1.6, 1.02), loc='upper right', title="Reactions")

    plt.tight_layout()

    plt.savefig(f'{dataset_name}.png', dpi=600)

    hyperdata_sorted.to_csv(f'{dataset_name}.csv', index=False)

def remove_exchanges(df):
    return df.drop([col for col in df.columns if col.startswith("EX_") or col.startswith("e_")], axis=1)

if __name__ == '__main__':
    os.chdir("../data/omics/")
    # achr_sample(r"light/HL/Dsalina_HL_Local2_2_4_4_fastcore_t0.xml", "e_Biomass__cytop")
    # achr_sample(r"light/ML/Dsalina_ML_Local2_2_4_4_fastcore_t0.xml", "e_Biomass__cytop")
    # achr_sample(r"light/LL/Dsalina_LL_Local2_2_4_4_fastcore_t0.xml", "e_Biomass__cytop")
    # achr_sample(r"nacl_h2o2_sorb/control/Dsalina_control_Local2_2_4_4_fastcore_t0.xml", "e_Biomass__cytop")
    # achr_sample(r"nacl_h2o2_sorb/nacl/Dsalina_nacl_Local2_2_4_4_fastcore_t0.xml", "e_Biomass__cytop")
    # achr_sample(r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_Local2_2_4_4_fastcore_t0.xml", "e_Biomass__cytop")
    # achr_sample(r"nacl_h2o2_sorb/sorb/Dsalina_sorb_Local2_2_4_4_fastcore_t0.xml", "e_Biomass__cytop")
    if not os.path.exists("pathways_map.csv"):
        model = MyModel(r"../models/model_with_trials.xml", "e_Biomass__cytop")
        model.get_pathway_reactions_map()
        df = pd.DataFrame.from_dict(data=model.pathway_reactions_map, orient='index').T
        df.to_csv("pathways_map.csv", index=False)
    # filenames = [r"light/HL_fastcore.xml", r"light/ML_fastcore.xml", r"light/LL_fastcore.xml"]
    # # Parallel(n_jobs=3)(delayed(achr_sample)(filename, "e_Biomass__cytop") for filename in filenames)
    # hl_samples = pd.read_csv(r"light/HL_fastcore_ACHR_samples.csv", index_col=0)
    # ml_samples = pd.read_csv(r"light/ML_fastcore_ACHR_samples.csv", index_col=0)
    # ll_samples = pd.read_csv(r"light/LL_fastcore_ACHR_samples.csv", index_col=0)
    # hl_samples.drop([col for col in hl_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    # ml_samples.drop([col for col in ml_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    # ll_samples.drop([col for col in ll_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    # hl_ml_results = kstest(hl_samples, ml_samples, "hl_ml")
    # ml_ll_results = kstest(ml_samples, ll_samples, "ml_ll")
    # pathway_enrichment(hl_ml_results, "hl_ml")
    # pathway_enrichment(ml_ll_results, "ml_ll")

    # control_samples = pd.read_csv(r"nacl_h2o2_sorb/control/Dsalina_control_Local2_2_4_4_fastcore_t0_ACHR_samples.csv", index_col=0)
    # nacl_samples = pd.read_csv(r"nacl_h2o2_sorb/nacl/Dsalina_nacl_Local2_2_4_4_fastcore_t0_ACHR_samples.csv", index_col=0)
    # sorb_samples = pd.read_csv(r"nacl_h2o2_sorb/sorb/Dsalina_sorb_Local2_2_4_4_fastcore_t0_ACHR_samples.csv", index_col=0)
    # h2o2_samples = pd.read_csv(r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_Local2_2_4_4_fastcore_t0_ACHR_samples.csv", index_col=0)
    # control_samples.drop([col for col in control_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    # nacl_samples.drop([col for col in nacl_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    # sorb_samples.drop([col for col in sorb_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    # h2o2_samples.drop([col for col in h2o2_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    # nacl_results = kstest(control_samples, nacl_samples, "nacl")
    # sorb_results = kstest(control_samples, sorb_samples, "sorb")
    # h2o2_results = kstest(control_samples, h2o2_samples, "h2o2")
    # pathway_enrichment(nacl_results, "nacl")
    # pathway_enrichment(sorb_results, "sorb")
    # pathway_enrichment(h2o2_results, "h2o2")

    filenames = [r"nacl_h2o2_sorb/control/Dsalina_control_gimme.xml", r"nacl_h2o2_sorb/nacl/Dsalina_nacl_gimme.xml", r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_gimme.xml", r"nacl_h2o2_sorb/sorb/Dsalina_sorb_gimme.xml" ]
    # achr_sample(r"nacl_h2o2_sorb/control/Dsalina_control_gimme.xml", "e_Biomass__cytop")
    # achr_sample(r"nacl_h2o2_sorb/nacl/Dsalina_nacl_gimme.xml", "e_Biomass__cytop")
    # achr_sample(r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_gimme.xml", "e_Biomass__cytop")
    # achr_sample(r"nacl_h2o2_sorb/sorb/Dsalina_sorb_gimme.xml", "e_Biomass__cytop")

    # Parallel(n_jobs=4)(delayed(achr_sample)(filename, "e_Biomass__cytop") for filename in filenames)

    control_samples = pd.read_csv(r"nacl_h2o2_sorb/control/Dsalina_control_gimme_ACHR_samples.csv", index_col=0)
    nacl_samples = pd.read_csv(r"nacl_h2o2_sorb/nacl/Dsalina_nacl_gimme_ACHR_samples.csv", index_col=0)
    sorb_samples = pd.read_csv(r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_gimme_ACHR_samples.csv", index_col=0)
    h2o2_samples = pd.read_csv(r"nacl_h2o2_sorb/sorb/Dsalina_sorb_gimme_ACHR_samples.csv", index_col=0)
    control_samples = remove_exchanges(control_samples)
    nacl_samples = remove_exchanges(nacl_samples)
    sorb_samples = remove_exchanges(sorb_samples)
    h2o2_samples = remove_exchanges(h2o2_samples)
    nacl_results = kstest(control_samples, nacl_samples, "nacl")
    sorb_results = kstest(control_samples, sorb_samples, "sorb")
    h2o2_results = kstest(control_samples, h2o2_samples, "h2o2")
    pathway_enrichment(nacl_results, "nacl")
    pathway_enrichment(sorb_results, "sorb")
    pathway_enrichment(h2o2_results, "h2o2")
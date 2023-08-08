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
from GSMMutils.model.COBRAmodel import MyModel
from GSMMutils import DATA_PATH
import logging
logging.getLogger('cobra.io').setLevel(logging.CRITICAL)


def split_reversible_reactions(model_to_sample):
    exchanges_demands_sinks = [reaction.id for reaction in model_to_sample.exchanges] + [reaction.id for reaction in model_to_sample.demands] + [reaction.id for reaction in model_to_sample.sinks]
    exchanges_demands_sinks = set(exchanges_demands_sinks)
    new_reactions = []
    for reaction in model_to_sample.reactions:
        if reaction not in exchanges_demands_sinks:
            if reaction.lower_bound < 0 < reaction.upper_bound:
                new_reaction = reaction.copy()
                new_reaction.id = reaction.id + "_reverse"
                new_reaction.lower_bound = 0
                new_reaction.upper_bound = -reaction.lower_bound
                for metabolite, coefficient in new_reaction.metabolites.items():
                    new_reaction.add_metabolites({metabolite: -coefficient})
                    new_reaction.add_metabolites({metabolite: -coefficient})
                new_reactions.append(new_reaction)
                reaction.lower_bound = 0
    model_to_sample.add_reactions(new_reactions)
    return model_to_sample


def achr_sample(filename, biomass_reaction):
    try:
        print('Sampling model: ', filename)
        folder = '/'.join(filename.split("/")[:-1])
        model_to_sample = read_sbml_model(filename)
        model_to_sample.objective = biomass_reaction
        initial_solution = model_to_sample.optimize().objective_value
        model_to_sample = split_reversible_reactions(model_to_sample)
        solution_after_spliting_reversible_reactions = model_to_sample.optimize().objective_value
        assert round(initial_solution, 5) == round(solution_after_spliting_reversible_reactions, 5)
        model_to_sample.reactions.get_by_id(biomass_reaction).lower_bound = 0.10
        for exchange in model_to_sample.exchanges:
            if exchange.lower_bound < 0 and not exchange.id.startswith("EX_C00205"):
                exchange.lower_bound = -10
        sampler = ACHRSampler(model_to_sample, thinning=100, seed=42)
        samples = sampler.sample(10000)
        result_filename = f"{folder}/{filename.split('/')[-1].split('.xml')[0]}_ACHR_samples.csv"
        samples.to_csv(result_filename, index=False)
    except Exception as e:
        print(e)
        print(f"Error in {filename}")


def load_results(file_names):
    data_list = []
    for filename in file_names:
        dataframe = pd.read_csv(filename)
        if "Unnamed: 0" in dataframe.columns:
            dataframe.drop(columns=["Unnamed: 0"], inplace=True)
        # dataframe = dataframe * dataframe[biomass_map[filename.split(".")[0]]]
        data_list.append(dataframe)
    reactions = list(set([reaction for dataframe in data_list for reaction in dataframe.columns]))
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


def kstest(samples_control: pd.DataFrame, samples_condition: pd.DataFrame, dataset_name: str):
    """
    Calculate the K-S test to detect significantly altered reactions fluxes.
    Results are saved in a csv file.

    Parameters
    ----------
    samples_control: pd.DataFrame
        The samples of the healthy tissue.
    samples_condition: pd.DataFrame
        The samples of the infected tissue.
    dataset_name: str
        The name of the dataset.
    Returns
    -------
    pd.DataFrame: The results of the K-S test for each reaction.

    """
    print("K-S test")
    union = set(samples_condition.columns).union(set(samples_control.columns))
    samples_condition_dict, samples_control_dict = {}, {}
    for rxn in union:
        if rxn not in samples_condition.columns:
            samples_condition_dict[rxn] = pd.Series(np.zeros(10000))
        if rxn not in samples_control.columns:
            samples_control_dict[rxn] = pd.Series(np.zeros(10000))

    samples_condition = pd.concat([samples_condition, pd.DataFrame.from_dict(samples_condition_dict)], axis=1)
    samples_control = pd.concat([samples_control, pd.DataFrame.from_dict(samples_control_dict)], axis=1)

    rxns1 = set(samples_condition.columns)
    rxns2 = set(samples_control.columns)

    rxn_c = rxns1.intersection(rxns2)
    # rxn_c = rxns1.union(rxns2)
    pvals = []
    rxnid = []
    fc = []

    for rxn in rxn_c:
        data1 = samples_condition[rxn].round(decimals=4)
        data2 = samples_control[rxn].round(decimals=4)

        data1 = data1.sample(n=1000)
        data2 = data2.sample(n=1000)

        if (data1.std() != 0 and data1.mean() != 0) or (data2.std() != 0 and data2.mean() != 0):
            kstat, pval = ks_2samp(data1, data2)

            data_1_mean = data1.mean()
            data_2_mean = data2.mean()

            foldc = (data_1_mean - data_2_mean) / abs(data_1_mean + data_2_mean)

            if data_1_mean < 0 and data_2_mean < 0:
                foldc = -foldc

            # if (data_1_mean < 0 < data_2_mean) or (data_1_mean > 0 > data_2_mean):
            #     print()

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

    rxns1 = set(samples_condition.columns)
    rxns2 = set(samples_control.columns)

    rxn_in1 = rxns1.difference(rxns2)
    rxn_in2 = rxns2.difference(rxns1)

    sigs = Parallel(n_jobs=8)(delayed(bootstrap_ci)(samples_condition[rx]) for rx in rxn_in1)
    act = [sigs[i][0] for i in range(len(sigs)) if sigs[i][1] == 1]
    sigs2 = Parallel(n_jobs=8)(delayed(bootstrap_ci)(samples_control[rx]) for rx in rxn_in2)
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
    try:
        print("Pathway enrichment...")
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
        to_ignore = ['Transporters pathway', 'Biosynthesis of cofactors', 'Microbial metabolism in diverse environments']  # 'Drains pathway',
        hyperdata_sorted = hyperdata_sorted[~hyperdata_sorted['Pathways'].isin(to_ignore)]

        plt.figure(figsize=(10, 5))

        sc = plt.scatter(hyperdata_sorted['P-value_adj'], np.arange(0, len(hyperdata_sorted['Pathways'])), s=hyperdata_sorted['ListReactions'], color=(0.9, 0.3, 0.1, 0.9))

        plt.xlabel('Adjusted p-value')

        plt.yticks(np.arange(0, len(hyperdata_sorted['Pathways'])), labels=hyperdata_sorted['Pathways'])

        handles, labels = sc.legend_elements(prop="sizes", alpha=0.8)

        plt.legend(handles, labels, bbox_to_anchor=(1.6, 1.02), loc='upper right', title="Reactions")

        plt.tight_layout()

        plt.savefig(f'{dataset_name}.png', dpi=600)

        hyperdata_sorted.to_csv(f'{dataset_name}.csv', index=False)
    except Exception as e:
        print(e)


def remove_exchanges(dataframe: pd.DataFrame) -> pd.DataFrame:
    return dataframe.drop([col for col in dataframe.columns if col.startswith("EX_") or col.startswith("e_")], axis=1)


if __name__ == '__main__':
    os.chdir(f"{DATA_PATH}/omics/")
    if not os.path.exists("pathways_map.csv"):
        model = MyModel(r"../models/model_with_trials.xml", "e_Biomass__cytop")
        model.get_pathway_reactions_map()
        results_dataframe = pd.DataFrame.from_dict(data=model.pathway_reactions_map, orient='index').T
        results_dataframe.to_csv("pathways_map.csv", index=False)
    filenames = [r"light/HL/Dsalina_HL_Local2_2_4_4_fastcore_t2.xml", r"light/ML/Dsalina_ML_Local2_2_4_4_fastcore_t2.xml", r"light/LL/Dsalina_LL_Local2_2_4_4_fastcore_t2.xml"]
    Parallel(n_jobs=3)(delayed(achr_sample)(filename, "e_Biomass__cytop") for filename in filenames)
    hl_samples = pd.read_csv(r"light/HL/Dsalina_HL_Local2_2_4_4_fastcore_t2_ACHR_samples.csv", index_col=0)
    ml_samples = pd.read_csv(r"light/ML/Dsalina_ML_Local2_2_4_4_fastcore_t2_ACHR_samples.csv", index_col=0)
    ll_samples = pd.read_csv(r"light/LL/Dsalina_LL_Local2_2_4_4_fastcore_t2_ACHR_samples.csv", index_col=0)
    hl_samples.drop([col for col in hl_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    ml_samples.drop([col for col in ml_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    ll_samples.drop([col for col in ll_samples.columns if col.startswith("EX_") or col.startswith("e_")], axis=1, inplace=True)
    hl_ml_results = kstest(hl_samples, ml_samples, "hl_ml")
    ml_ll_results = kstest(ml_samples, ll_samples, "ml_ll")
    hl_ll_results = kstest(hl_samples, ll_samples, "hl_ll")
    pathway_enrichment(hl_ml_results, "hl_ml")
    pathway_enrichment(ml_ll_results, "ml_ll")
    pathway_enrichment(hl_ll_results, "hl_ll")

    # filenames = [
        # r"nacl_h2o2_sorb/control/Dsalina_control_gimme.xml",
        # r"nacl_h2o2_sorb/nacl/Dsalina_nacl_gimme.xml",
        # r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_gimme.xml",
        # r"nacl_h2o2_sorb/sorb/Dsalina_sorb_gimme.xml",
    #     r"nacl_h2o2_sorb/control/Dsalina_control_Local2_2_4_4_fastcore_t3_7.xml",  # nacl_h2o2_sorb/control/Dsalina_control_Local2_2_4_4_fastcore_t3_7.xml
    #     r"nacl_h2o2_sorb/nacl/Dsalina_nacl_Local2_2_4_4_fastcore_t2_8.xml",
    #     r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_Local2_2_4_4_fastcore_t3_7.xml",
    #     r"nacl_h2o2_sorb/sorb/Dsalina_sorb_Local2_2_4_4_fastcore_t2_8.xml"
    # ]
    # Parallel(n_jobs=8)(delayed(achr_sample)(filename, "e_Biomass__cytop") for filename in filenames)
    #
    # control_samples_gimme = pd.read_csv(r"nacl_h2o2_sorb/control/Dsalina_control_gimme_ACHR_samples.csv")
    # nacl_samples_gimme = pd.read_csv(r"nacl_h2o2_sorb/nacl/Dsalina_nacl_gimme_ACHR_samples.csv")
    # sorb_samples_gimme = pd.read_csv(r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_gimme_ACHR_samples.csv")
    # h2o2_samples_gimme = pd.read_csv(r"nacl_h2o2_sorb/sorb/Dsalina_sorb_gimme_ACHR_samples.csv")
    #
    # control_samples_gimme = remove_exchanges(control_samples_gimme)
    # nacl_samples_gimme = remove_exchanges(nacl_samples_gimme)
    # sorb_samples_gimme = remove_exchanges(sorb_samples_gimme)
    # h2o2_samples_gimme = remove_exchanges(h2o2_samples_gimme)
    #
    # nacl_results_gimme = kstest(control_samples_gimme, nacl_samples_gimme, "nacl_gimme")
    # sorb_results_gimme = kstest(control_samples_gimme, sorb_samples_gimme, "sorb_gimme")
    # h2o2_results_gimme = kstest(control_samples_gimme, h2o2_samples_gimme, "h2o2_gimme")
    # pathway_enrichment(nacl_results_gimme, "nacl_gimme")
    # pathway_enrichment(sorb_results_gimme, "sorb_gimme")
    # pathway_enrichment(h2o2_results_gimme, "h2o2_gimme")
    #
    # control_samples_fastcore = pd.read_csv(r"nacl_h2o2_sorb/control/Dsalina_control_Local2_2_4_4_fastcore_t3_7_ACHR_samples.csv", index_col=0)
    # nacl_samples_fastcore = pd.read_csv(r"nacl_h2o2_sorb/nacl/Dsalina_nacl_Local2_2_4_4_fastcore_t2_8_ACHR_samples.csv", index_col=0)
    # sorb_samples_fastcore = pd.read_csv(r"nacl_h2o2_sorb/h2o2/Dsalina_h2o2_Local2_2_4_4_fastcore_t3_7_ACHR_samples.csv", index_col=0)
    # h2o2_samples_fastcore = pd.read_csv(r"nacl_h2o2_sorb/sorb/Dsalina_sorb_Local2_2_4_4_fastcore_t2_8_ACHR_samples.csv", index_col=0)
    # control_samples_fastcore = remove_exchanges(control_samples_fastcore)
    # nacl_samples_fastcore = remove_exchanges(nacl_samples_fastcore)
    # sorb_samples_fastcore = remove_exchanges(sorb_samples_fastcore)
    # h2o2_samples_fastcore = remove_exchanges(h2o2_samples_fastcore)
    #
    # nacl_results_fastcore = kstest(control_samples_fastcore, nacl_samples_fastcore, "nacl_fastcore")
    # sorb_results_fastcore = kstest(control_samples_fastcore, sorb_samples_fastcore, "sorb_fastcore")
    # h2o2_results_fastcore = kstest(control_samples_fastcore, h2o2_samples_fastcore, "h2o2_fastcore")
    # pathway_enrichment(nacl_results_fastcore, "nacl_fastcore")
    # pathway_enrichment(sorb_results_fastcore, "sorb_fastcore")
    # pathway_enrichment(h2o2_results_fastcore, "h2o2_fastcore")

    print("done")

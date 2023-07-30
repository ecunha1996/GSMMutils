import os
from os.path import join
from typing import Union

import matlab.engine
import numpy as np
import pandas as pd
import scipy
from cobra import Model

from GSMMutils import SRC_PATH, DATA_PATH
from GSMMutils.model.COBRAmodel import MyModel


def normalize_to_one(dataframe):
    dataframe['value'] = dataframe['value'] / dataframe['value'].sum()
    return dataframe.round(5)


def parse_other_lipids(lipid_abb, list_of_names, chains_map):
    final_map = {}
    for lipid in list_of_names:
        try:
            acyl = lipid.split("-sn-")[0].split("-3-")[0].split("-(6'-")[0]
            acyl = acyl.replace('5Z,8Z,11Z,14Z-', "").replace('8Z,11Z,14Z-', "").replace('8Z,11Z,14Z,17Z-',
                                                                                         "").replace('11Z,14Z,17Z-',
                                                                                                     "").replace(
                "-3-beta-D-galactosyl", "").replace("-3-O-?-D-galactosyl", "").replace("7Z,10Z-", "").replace(
                "-3-O-D-galactosyl", "")
            if "1-" in acyl and "-2-" in acyl:
                chain_1 = acyl.split("-2-")[0].lstrip("1-").strip("(").strip(")")
                chain_2 = acyl.split("-2-")[1].strip("(").strip(")")
                if "6Z" not in chain_1:
                    chain_1 = chain_1.replace("9Z,12Z-", "")
                if "6Z" not in chain_2:
                    chain_2 = chain_2.replace("9Z,12Z-", "")
                final_map[lipid] = lipid_abb + "__" + chains_map[chain_1] + "__" + chains_map[chain_2]
            elif "-di" in acyl:
                chains = acyl.split("-di")[1].lstrip("-").strip("(").strip(")")
                if "6Z" not in chains:
                    chains = chains.replace("9Z,12Z-", "")
                chain_1 = chains
                chain_2 = chains
                final_map[lipid] = lipid_abb + "__" + chains_map[chain_1] + "__" + chains_map[chain_2]
            else:
                print(acyl)
        except Exception as e:
            print(e)
    as_list = list(final_map.values())
    for e in as_list:
        if as_list.count(e) > 1:
            print(e)
    return final_map


def get_metmat_tag(final_map, set_of_names, chains_map):
    faas = []
    met_mat = [[1 for _ in range(len(final_map.keys()))]]
    for fatty_acid_name, code in chains_map.items():
        temp = []
        for lipid in set_of_names:
            splited = lipid.split("__")
            chain_1 = splited[1]
            chain_2 = splited[2]
            chain_3 = splited[3]
            counter = [chain_1, chain_2, chain_3]
            temp.append(counter.count(code))
            faas += counter
        met_mat.append(temp)
    return met_mat, faas


def get_metmat(final_map, set_of_names, chains_map):
    faas = []
    met_mat = [[1 for _ in range(len(final_map.keys()))]]
    for fatty_acid_name, code in chains_map.items():
        temp = []
        for lipid in set_of_names:
            splited = lipid.split("__")
            chain_1 = splited[1]
            chain_2 = splited[2]
            counter = [chain_1, chain_2]
            temp.append(counter.count(code))
            faas += counter
        met_mat.append(temp)
    return met_mat, faas


def parse_dgts(model):
    list_of_names = [e.name for e in model.reactions.e_DGTS_complete__e_r_.reactants]
    final_map = {}
    for lipid in list_of_names:
        acyl = lipid
        chain_1 = acyl.split("/")[0].split(" ")[1].lstrip("(").replace(":", "_")
        if "6Z" not in chain_1:
            chain_1 = chain_1.split("(")[0]
        else:
            chain_1 = chain_1.split("(")[0] + "v2"
        chain_2 = acyl.split("/")[1].rstrip(")").replace(":", "_")
        if "6Z" not in chain_2:
            chain_2 = chain_2.split("(")[0]
        else:
            chain_2 = chain_2.split("(")[0] + "v2"
        final_map[lipid] = "DGTS__" + chain_1 + "__" + chain_2
    return final_map


def parse_tag(model):
    list_of_names = [e.name for e in model.reactions.e_TAG_complete__lip.reactants]
    chains_map = {'dodecanoyl': "12_0", "tetradecanoyl": "14_0", '9Z-tetradecenoyl': "14_1", "hexadecanoyl": "16_0",
                  "9Z-hexadecenoyl": "16_1", '9Z,12Z-hexadecadienoyl': '16_2',
                  "9Z,12Z,15Z-hexadecatrienoyl": "16_3", "4Z,7Z,10Z,13Z-hexadecatetraenoyl": "16_4",
                  "octadecanoyl": "18_0", "9Z-octadecenoyl": "18_1", '9Z,12Z-octadecadienoyl': "18_2",
                  "6Z,9Z,12Z-octadecatrienoyl": '18_3v2', '9Z,12Z,15Z-octadecatrienoyl': "18_3",
                  'eicosanoyl': "20_0", '11Z-eicosenoyl': "20_1", '11Z,14Z-eicosadienoyl': "20_2",
                  'eicosatrienoyl': "20_3", 'eicosatetraenoyl': "20_4",
                  "5Z,8Z,11Z,14Z,17Z-eicosapentaenoyl": "20_5", "docosanoyl": "22_0", "13Z-docosenoyl": "22_1",
                  '13Z,16Z-docosadienoyl': "22_2", "4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl": "22_6"}
    final_map = {}
    for lipid in list_of_names:
        try:
            if "eicosapentaenoyl" in lipid and "5Z,8Z,11Z,14Z,17Z-eicosapentaenoyl" not in lipid:
                print()
            acyl = lipid.split("-sn-")[0]
            acyl = acyl.replace('(8Z,11Z,14Z-', "").replace('(5Z,8Z,11Z,14Z-', "").replace('(8Z,11Z,14Z,17Z-',
                                                                                           "").replace(
                '(11Z,14Z,17Z', "")
            if "1-" in acyl and "-2-" in acyl and "-3-" in acyl:
                chain_1 = acyl.split("-2-")[0].lstrip("1-").strip(")").strip("(")
                chain_2 = acyl.split("-2-")[1].split("-3-")[0].strip(")").strip("(")
                chain_3 = acyl.split("-3-")[1].strip(")").strip("(")
                final_map[lipid] = "TAG__" + chains_map[chain_1] + "__" + chains_map[chain_2] + "__" + chains_map[
                    chain_3]
            elif "1,2-di" in acyl:
                chains = acyl.split("-3-")[0].replace("1,2-di", "").lstrip("-").strip(")").strip("(")
                chain_3 = acyl.split("-3-")[1].strip(")").strip("(")
                chain_1 = chains
                chain_2 = chains
                final_map[lipid] = "TAG__" + chains_map[chain_1] + "__" + chains_map[chain_2] + "__" + chains_map[
                    chain_3]
            elif "2,3-di" in acyl:
                chain_1 = acyl.split("-2")[0].replace("1-", "").strip(")").strip("(")
                chains = acyl.split("3-")[1].strip("di").lstrip("-").strip(")").strip("(")
                chain_2 = chains
                chain_3 = chains
                final_map[lipid] = "TAG__" + chains_map[chain_1] + "__" + chains_map[chain_2] + "__" + chains_map[
                    chain_3]
            elif "1,3-di" in acyl:
                chains = acyl.split("-2-")[0].strip("1,3-di").lstrip("-").strip(")").strip("(")
                chain_1 = chains
                chain_3 = chains
                chain_2 = acyl.split("-2-")[1].strip(")").strip("(")
                final_map[lipid] = "TAG__" + chains_map[chain_1] + "__" + chains_map[chain_2] + "__" + chains_map[
                    chain_3]
            elif "1,2,3-tri" in acyl:
                chains = acyl.replace("1,2,3-tri", "").strip("-").strip(")").strip("(")
                chain_1 = chains
                chain_2 = chains
                chain_3 = chains
                final_map[lipid] = "TAG__" + chains_map[chain_1] + "__" + chains_map[chain_2] + "__" + chains_map[
                    chain_3]
            else:
                print(acyl)
        except Exception as e:
            print(e)
            print(lipid)
            print("----")
    as_list = list(final_map.values())
    for e in as_list:
        if as_list.count(e) > 1:
            print(e)
    return final_map


def parse_2fa_lipid(lipid_abb, model, compartment_id="C_00003", parent_reaction=None):
    compartment = model.compartments[compartment_id]
    if not parent_reaction:
        parent_reaction = f'e_{lipid_abb}__{compartment}'
    list_of_names = [e.name for e in model.reactions.get_by_id(parent_reaction).reactants]
    chains_map = {'dodecanoyl': "12_0", "tetradecanoyl": "14_0", '9Z-tetradecenoyl': "14_1", "hexadecanoyl": "16_0",
                  "9Z-hexadecenoyl": "16_1", 'hexadecadienoyl': '16_2', "9Z,12Z,15Z-hexadecatrienoyl": "16_3",
                  "4Z,7Z,10Z,13Z-hexadecatetraenoyl": "16_4", "octadecanoyl": "18_0", "9Z-octadecenoyl": "18_1",
                  'octadecadienoyl': "18_2", "6Z,9Z,12Z-octadecatrienoyl": '18_3v2',
                  '9Z,12Z,15Z-octadecatrienoyl': "18_3", 'eicosanoyl': "20_0", '11Z-eicosenoyl': "20_1",
                  '11Z,14Z-eicosadienoyl': "20_2", 'eicosatrienoyl': "20_3", 'eicosatetraenoyl': "20_4",
                  "eicosapentaenoyl": "20_5", "docosanoyl": "22_0", "13Z-docosenoyl": "22_1",
                  '13Z,16Z-docosadienoyl': "22_2", "4Z,7Z,10Z,13Z,16Z,19Z-docosahexaenoyl": "22_6"}
    if lipid_abb == "DGTS":
        final_map = parse_dgts(model)
    elif lipid_abb == "TAG":
        final_map = parse_tag(model)
    else:
        final_map = parse_other_lipids(lipid_abb, list_of_names, chains_map)

    set_of_names = set(final_map.values())
    if lipid_abb == "TAG":
        met_mat, faas = get_metmat_tag(final_map, set_of_names, chains_map)
    else:
        met_mat, faas = get_metmat(final_map, set_of_names, chains_map)
    if lipid_abb == "DAG":
        set_of_names = list(set_of_names) + ['DAG__18_3__18_1__chlo', 'DAG__18_3v2__16_0__chlo',
                                             'DAG__18_3__16_0__chlo', 'DAG__18_2__18_1__chlo',
                                             'DAG__18_2__16_0__chlo', 'DAG__18_1__18_1__chlo',
                                             'DAG__18_1__16_0__chlo']
    # print(f'lipidClass = {set_of_names};')
    i = 0
    chlo_matrix = [[1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0],
                   [0, 1, 1, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0],
                   [1, 0, 0, 1, 0, 2, 1], [0, 0, 0, 1, 1, 0, 0], [0, 1, 0, 0, 0, 0, 0], [1, 0, 1, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0]]
    final_met_mat = []
    for row in met_mat:
        if not all(e == 0 for e in row):
            if lipid_abb == "DAG":
                temp_row = row + chlo_matrix[i]
            else:
                temp_row = row
            final_met_mat.append(temp_row)
            # print(','.join (str(e) for e in temp_row) + ";")
            i += 1
    final_map_rev = {}
    for key, value in final_map.items():
        final_map_rev[value] = model.get_metabolite_by_name(key, compartment=compartment_id).id.split("__")[0]
    return final_met_mat, final_map, faas, final_map_rev, set_of_names


def merge_results():
    directory = r"C:\Users\Bisbii\OneDrive - Universidade do Minho\Algae\Models\Dsalina\lipids"
    os.chdir(directory)
    res = pd.DataFrame()
    for file in os.listdir("/"):
        if file.endswith(".csv") and file != "model_stoichiometry_all.csv":
            results_dataframe = pd.read_csv(file, sep="\t")
            res = pd.concat([res, results_dataframe], axis=0)
    print(res.shape)
    res = res.drop_duplicates()
    print(res.shape)
    res.to_csv("model_stoichiometry_all.csv", index=False, sep="\t")


class FAME2Biomass:
    def __init__(self, model: Union[MyModel, Model], data_directory):
        self.model = model
        self.data_directory = data_directory

    def parse_results_to_merlin(self, final_map_rev, reaction_id, lipid='pg'):
        res = []
        compounds = pd.read_csv(rf"{self.data_directory}/model_ngaditana_lipids/model_compound_from_db.csv")
        lipid_results = pd.read_csv(
            f"{self.data_directory}/model_ngaditana_lipids/{lipid}_opt_results.tsv",
            sep="\t", index_col=0,
            header=None).astype(float)
        lipid_results.columns = ["value"]
        lipid_results = normalize_to_one(lipid_results)
        lipid_results.to_csv(f"{lipid}_normalized.tsv", sep="\t")
        as_dict = lipid_results.to_dict(orient='index')
        compartment = 14 if lipid.upper() != "TAG" else 57
        if lipid.upper() in ("PG", "MGDG", "DGDG", "SQDG"):
            compartment = 13
        for key, value in as_dict.items():
            if "chlo" in key:
                key = key.replace("__chlo", "")
                compartment = 13
            idcompound = compounds.loc[compounds["external_identifier"] == final_map_rev[key], 'idcompound'].values[0]
            res.append([-value["value"], compartment, idcompound, reaction_id])
        res_df = pd.DataFrame(res, columns=["stoichiometry", "compartment_id", "idcompound", "idreaction"])
        res_df.to_csv(f"{self.data_directory}/model_ngaditana_lipids/{lipid}_model_stoichiometry.csv", index=False)

    def run(self, lipid_abb, compartment="er", compartment_id="C_00001"):
        met_mat, _, _, final_map_rev, lipid_class = parse_2fa_lipid(lipid_abb, self.model,
                                                                    compartment_id=compartment_id,
                                                                    parent_reaction=f"e_{lipid_abb}_complete"
                                                                                    f"__{compartment}")
        eng = matlab.engine.start_matlab()
        eng.cd(self.data_directory)
        eng.addpath(join(SRC_PATH, 'matlab'))
        eng.addpath(r'C:\gurobi912\win64\matlab')
        b_list = [1 / 3, 0, 0.066765796, 0, 0.393942655, 0.415415317, 0.018880291, 0.004891267, 0, 0.019809632,
                  0.042309461, 0.007630377, 0.004891267 / 2, 0.004891267 / 2, 0, 0.003326062, 0.000821733, 0.000821733,
                  0.008999932, 0.011494478, 0, 0, 0, 0]
        b_vector = matlab.double(vector=[e * 3 for e in b_list])
        np_a = np.array(met_mat)
        scipy.io.savemat(rf"{DATA_PATH}\fame2biomass\{self.model.id}\{lipid_abb.lower()}_metmat.mat", {'metmat': np_a})
        function_name = getattr(eng, f"lipids_{lipid_abb.lower()}")
        function_name(list(lipid_class), b_vector)
        # eng.fame2model(list(lipid_class), lipid_abb.lower())
        os.chdir(self.data_directory)
        reactions_ids_map = {"PG": 35724, "PE": 36207, "PI": 35940, "PC": 35921, "MGDG": 35956, "DGDG": 35957,
                             "SQDG": 35962, "DGTS": 35955, "DAG": 36091, "TAG": 52720}
        self.parse_results_to_merlin(final_map_rev, reactions_ids_map[lipid_abb], lipid=lipid_abb.lower())

    def batch_run(self, lipid_compartments_map, abb_compartment_id_map):
        for lipid, compartment in lipid_compartments_map.items():
            self.run(lipid, compartment, abb_compartment_id_map[compartment])

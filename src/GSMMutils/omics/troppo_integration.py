import json
import sys
import os
import warnings
import traceback
from os.path import abspath, dirname, join

import pandas
import pandas as pd
from numpy import linspace

from GSMMutils import DATA_PATH

# Warnings and prints are blocked to avoid massive outputs that appear after running COBAMP.
warnings.simplefilter("ignore")

CONFIG_PATH = abspath(join(dirname(__file__), '../../../config'))
# Disable
def block_print():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enable_print():
    sys.stdout = sys.__stdout__


block_print()

import cobra
import re

from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from cobamp.utilities.parallel import batch_run

enable_print()
#
# import time
# while True:
#     time.sleep(10)

params = json.load(open(rf"{CONFIG_PATH}/troppo_nacl.json", "r"))
media = pd.read_excel(rf"{DATA_PATH}/media.xlsx", index_col=0, sheet_name=None, engine='openpyxl')['media_with_starch'].to_dict(orient='index')
MEDIUM_CONDITIONS = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}

UPTAKE_DRAINS= {'f2 medium': {"EX_C00014__dra", "EX_C00059__dra", "EX_C01330__dra", "EX_C00011__dra", "EX_C14818__dra", "EX_C00080__dra", "EX_C00001__dra", "EX_C00305__dra", "EX_C01327__dra", "EX_C00244__dra", "EX_C00009__dra",
                               "EX_C00007__dra", "EX_C00205__dra", "EX_C00378__dra", "EX_C00120__dra", "EX_C02823__dra", 'DM_C00369__chlo'
                               }}

MEDIUM_METABOLITES = {'f2 medium': [e.split("__")[0].replace("EX_", "") for e in MEDIUM_CONDITIONS['f2 medium']]}

#THRESHOLDS = linspace(0, 9, 10)  # List of ints or floats.
THRESHOLDS = [2.8, 3.7]
INTEGRATION_THRESHOLDS = list(THRESHOLDS)

AND_OR_FUNCS = (min, sum)

PROTECTED = ["e_Biomass__cytop",] + list(UPTAKE_DRAINS['f2 medium'])

MODEL_PATH = rf"{DATA_PATH}/models/model_with_trials.xml"
CONSISTENT_MODEL_PATH = rf"{DATA_PATH}/models/consistent_model.xml"
OMICS_DATA_PATH = rf"{DATA_PATH}/omics/raw_counts.txt"
TROPPO_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
MODEL_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
MODEL_TASKS_PATH = f"{DATA_PATH}/omics"
MODEL_TASKS_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"

def reconstruction_function(omics_container, parameters: dict):
    """
    This function is used to run the reconstruction algorithm.

    Parameters
    ----------
    omics_container : pandas.DataFrame
        The omics data set.
    parameters : dict
        The parameters to be used for the reconstruction algorithm.

    Returns
    ----------
    rec_wrapper : Reconstruction Wrapper object with model and omics data.
    """
    block_print()

    def integration_fx(data_map):
        return [[k for k, v in data_map.get_scores().items() if (v is not None and v > threshold) or k in PROTECTED]]

    def score_apply(data_map):
        dm = {k: 0 if v is None else (min(v, 10) - threshold) if k not in PROTECTED else 20
              for k, v in data_map.items()}
        return dm

    threshold, rec_wrapper, method = [parameters[parameter] for parameter in
                                      ['threshold', 'reconstruction_wrapper', 'algorithm']]

    # noinspection PyBroadException
    try:
        if method == 'fastcore':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=AND_OR_FUNCS,
                                              integration_strategy=('custom', [integration_fx]), solver='CPLEX')

        elif method == 'tinit':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=AND_OR_FUNCS,
                                              integration_strategy=('continuous', score_apply), solver='CPLEX')

    except Exception as e:
        traceback.print_exc()
        return {r: False for r in rec_wrapper.model_reader.r_ids}


def troppo_omics_integration(model: cobra.Model, algorithm: str, threshold: float, thread_number: int,
                             omics_dataset: pandas.DataFrame, thresholds_map=None):
    """
    This function is used to run the Troppo's integration algorithms.

    Parameters
    ----------
    omics_dataset: pandas.DataFrame
        A dataframe containing the omics dataset to integrate
    model: cobra.Model
        The COBRA model.
    algorithm: str
        The algorithm to be used.
    threshold: float
        The threshold to be used.
    thread_number: int
        The number of threads to be used.

    Returns
    -------
    integration_results: dict
        Dataframe containing the results of the omics integration.
        Each sample as a dictionary containing a boolean value for each reaction.
    """
    if thresholds_map is None:
        thresholds_map = {}
    block_print()

    patt = re.compile('__COBAMPGPRDOT__[0-9]{1}')
    replace_alt_transcripts = lambda x: patt.sub('', x)

    details = [params['DATASET'], algorithm, threshold]

    template = model.copy()

    omics_data = TabularReader(path_or_df=omics_dataset, nomenclature=params['NOMENCLATURE'],
                               omics_type=params['OMICS_TYPE']).to_containers()

    enable_print()
    print('Tabular Reader Finished.')
    block_print()

    reconstruction_wrapper = ReconstructionWrapper(model=template, ttg_ratio=9999,
                                                   gpr_gene_parse_function=replace_alt_transcripts)

    enable_print()
    print('Reconstruction Wrapper Finished.')
    # block_print()

    parameters = {'threshold': threshold, 'reconstruction_wrapper': reconstruction_wrapper, 'algorithm': algorithm}

    batch_fastcore_res = batch_run(reconstruction_function, omics_data, parameters, threads=thread_number)

    result_dict = dict(zip([sample.condition for sample in omics_data], batch_fastcore_res))

    enable_print()
    print(f'Omics Integration with {details[1]} (Threshold = {details[2]}) Finished.')

    return result_dict

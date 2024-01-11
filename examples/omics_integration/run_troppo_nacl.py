import json
import math
import os
import sys
from os.path import join, dirname, abspath

import pandas as pd
from numpy import linspace

from gsmmutils import DATA_PATH
from gsmmutils.omics.omics_processing import thresholding_filter
from gsmmutils.omics.troppo_integration import troppo_omics_integration
sys.path.insert(0, r"/")
sys.path.insert(0, "/home/src/")
from gsmmutils.io import read_csv
from gsmmutils.omics.omics_integration import OmicsIntegration
from gsmmutils.omics.model_handle import load_model, sbml_model_reconstruction
CONFIG_PATH = abspath(join(dirname(__file__), '../../config'))
params = json.load(open(rf"{CONFIG_PATH}/troppo_nacl.json", "r"))
MODEL_PATH = rf"{DATA_PATH}/models/model_with_trials.xml"
CONSISTENT_MODEL_PATH = rf"{DATA_PATH}/models/consistent_model.xml"
OMICS_DATA_PATH = rf"{DATA_PATH}/omics/raw_counts.txt"
TROPPO_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
MODEL_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
MODEL_TASKS_PATH = f"{DATA_PATH}/omics"
MODEL_TASKS_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
media = pd.read_excel(rf"{DATA_PATH}/media.xlsx", index_col=0, sheet_name=None, engine='openpyxl')['media_with_starch'].to_dict(orient='index')
MEDIUM_CONDITIONS = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}

MEDIUM_METABOLITES = {'f2 medium': [e.split("__")[0].replace("EX_", "") for e in MEDIUM_CONDITIONS['f2 medium']]}

#THRESHOLDS = linspace(0, 9, 10)  # List of ints or floats.
THRESHOLDS = [2.8, 3.7]
INTEGRATION_THRESHOLDS = list(THRESHOLDS)


def reconstruction_pipeline():
    """
    This function is used to run the reconstruction pipeline.

    In the end this function generates a SBML file of the tissue-specific reconstructed_models that resulted from the
    omics data integration with Troppo.

    All the parameters required to run this pipeline can be defined in the ***pipeline_paths.py***,
    ***config_variables.py***, and ***medium_variables.py files***.
    """
    print('-------------------------------------------------------------------------------------------------------')
    print('--------------------------------------- Loading Template Model. ---------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    template_model = load_model(model_path=MODEL_PATH, consistent_model_path=CONSISTENT_MODEL_PATH)

    print('-------------------------------------------------------------------------------------------------------')
    print('-------------------------------------- Processing Omics Dataset. --------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')
    omics = OmicsIntegration('omics/raw_counts.txt', samples_names={"SRR6825159_Aligned.sortedByCoord.out.bam": "control_1",
                                                                    "SRR6825160_Aligned.sortedByCoord.out.bam": "control_2",
                                                                    "SRR6825161_Aligned.sortedByCoord.out.bam": "control_3",
                                                                    "SRR6825162_Aligned.sortedByCoord.out.bam": "nacl_1",
                                                                    "SRR6825163_Aligned.sortedByCoord.out.bam": "nacl_2",
                                                                    "SRR6825164_Aligned.sortedByCoord.out.bam": "nacl_3",
                                                                    "SRR6825165_Aligned.sortedByCoord.out.bam": "h2o2_1",
                                                                    "SRR6825166_Aligned.sortedByCoord.out.bam": "h2o2_2",
                                                                    "SRR6825167_Aligned.sortedByCoord.out.bam": "h2o2_3",
                                                                    "SRR6825168_Aligned.sortedByCoord.out.bam": "sorb_1",
                                                                    "SRR6825169_Aligned.sortedByCoord.out.bam": "sorb_2",
                                                                    "SRR6825170_Aligned.sortedByCoord.out.bam": "sorb_3",
                                                                    }, model=template_model)
    omics.getmm = read_csv(r"omics/getmm_salinity.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics.sum_tech_reps()
    omics_data = omics.counts.applymap(lambda x: math.log2(x + 1))
    omics_data_temp = omics_data.loc[(omics_data != 0).any(axis=1)]
    print(omics_data_temp.describe())
    # import matplotlib.pyplot as plt
    # fig = omics_data.plot.density()
    # fig.set_xticks([x / 2 for x in range(0, 40)])
    # plt.savefig(r"omics/counts_density_nacl.png")
    omics_data = omics_data.T
    print('Omics dataset Loaded.')

    if params['THRESHOLDING_STRATEGY'] != 'default':
        omics_data = thresholding_filter(omics_dataframe=omics_data, thresholding_strategy=params['THRESHOLDING_STRATEGY'],
                                         global_threshold_upper=params['GLOBAL_THRESHOLD_UPPER'],
                                         global_threshold_lower=params['GLOBAL_THRESHOLD_LOWER'],
                                         local_threshold=params['LOCAL_THRESHOLD'])

        print(f'{params["THRESHOLDING_STRATEGY"]} threshold filter applied.')

    print('-------------------------------------------------------------------------------------------------------')
    print('------------------------------- Starting Omics Integration with Troppo. -------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    integration_result = {}

    for algorithm in params["ALGORITHMS"]:
        troppo_result_dict = {}
        print(f'Omics integration with {algorithm} started.')
        print('-------------------------------------------------------------------------------------------------------')

        if algorithm == 'fastcore':
            thred_number = params["THREAD_NUMBER_FASTCORE"]

        elif algorithm == 'tinit':
            thred_number = params["THREAD_NUMBER_TINIT"]

        else:
            return 'Algorithm not supported by this pipeline.'

        for threshold in INTEGRATION_THRESHOLDS:
            troppo_result = troppo_omics_integration(model=template_model, algorithm=algorithm, threshold=threshold,
                                                     thread_number=thred_number, omics_dataset=omics_data)

            for sample in list(troppo_result.keys()):
                th = str(round(threshold, 2)).replace('.', '_')
                troppo_result_dict[f'{params["MODEL"]}_{sample}_{algorithm}_t{th}'] = troppo_result[sample]
                integration_result[f'{params["MODEL"]}_{sample}_{algorithm}_t{th}'] = troppo_result[sample]

            print('----------------------------------------------------------'
                  '---------------------------------------------')

        if params["THRESHOLDING_STRATEGY"] == 'default':
            result_path = os.path.join(TROPPO_RESULTS_PATH,
                                       f'{params["MODEL"]}_{params["DATASET"]}_{algorithm}_{params["THRESHOLDING_STRATEGY"]}.csv')
        else:
            file_path = f'{params["MODEL"]}_{params["DATASET"]}_{algorithm}_{params["THRESHOLDING_STRATEGY"]}_' \
                        f'{params["GLOBAL_THRESHOLD_UPPER"]}_{params["GLOBAL_THRESHOLD_LOWER"]}_{params["LOCAL_THRESHOLD"]}.csv'
            result_path = os.path.join(TROPPO_RESULTS_PATH, file_path)

        troppo_result_dataframe = pd.DataFrame.from_dict(troppo_result_dict, orient='index')
        troppo_result_dataframe.to_csv(result_path)

    if params["EVALUATE_TASKS"]:
        print('------------------------ Starting task evaluation for the integration results. ------------------------')
        print('-------------------------------------------------------------------------------------------------------')

        # task_eval(model_template=template_model, integration_results=integration_result)

        print('-------------------------------------------------------------------------------------------------------')

    if params["RECONSTRUCT_MODELS"]:
        print('----------------------- Starting reconstruction of the context-specific models. -----------------------')
        print('-------------------------------------------------------------------------------------------------------')

        for sample_name in list(integration_result.keys()):
            print(f'Context-specific model reconstruction for {sample_name} started.')
            print('---------------------------------------------------------------'
                  '----------------------------------------')
            sbml_model_reconstruction(original_model=template_model, sample=sample_name,
                                      integration_result_dict=integration_result)

            print('---------------------------------------------------------------'
                  '----------------------------------------')

    print('------------------------------------------ Pipeline Finished ------------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')


if __name__ == '__main__':
    os.chdir(DATA_PATH)
    reconstruction_pipeline()

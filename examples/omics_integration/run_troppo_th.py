import json
import math
import os
from os.path import abspath, dirname, join
import pandas as pd
from GSMMutils import DATA_PATH
from GSMMutils.io import read_csv
from GSMMutils.omics.omics_integration import OmicsIntegration
from GSMMutils.omics.omics_processing import thresholding_filter
from GSMMutils.omics.troppo_integration import troppo_omics_integration
from GSMMutils.omics.model_handle import load_model, sbml_model_reconstruction

CONFIG_PATH = abspath(join(dirname(__file__), '../../config'))

params = json.load(open(rf"{CONFIG_PATH}/troppo_light.json", "r"))
UPTAKE_DRAINS = {'f2 medium': {"EX_C00014__dra", "EX_C00059__dra", "EX_C01330__dra", "EX_C00011__dra", "EX_C14818__dra", "EX_C00080__dra", "EX_C00001__dra", "EX_C00305__dra", "EX_C01327__dra", "EX_C00244__dra", "EX_C00009__dra",
                               "EX_C00007__dra", "EX_C00205__dra", "EX_C00378__dra", "EX_C00120__dra", "EX_C02823__dra", 'DM_C00369__chlo'
                               }}
params['uptake_drains'] = UPTAKE_DRAINS
MODEL_PATH = rf"{DATA_PATH}/models/model_with_trials.xml"
CONSISTENT_MODEL_PATH = rf"{DATA_PATH}/models/consistent_model_v2.xml"
OMICS_DATA_PATH = rf"{DATA_PATH}/omics/raw_data_light.txt"
TROPPO_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
MODEL_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
MODEL_TASKS_PATH = f"{DATA_PATH}/omics"
MODEL_TASKS_RESULTS_PATH = f"{DATA_PATH}/omics/{params['DATASET']}"
media = pd.read_excel(rf"{DATA_PATH}/media.xlsx", index_col=0, sheet_name=None, engine='openpyxl')['media_with_starch'].to_dict(orient='index')
MEDIUM_CONDITIONS = {'f2 medium': {key: (value['LB'], value['UB']) for key, value in media.items()}}

MEDIUM_METABOLITES = {'f2 medium': [e.split("__")[0].replace("EX_", "") for e in MEDIUM_CONDITIONS['f2 medium']]}

# THRESHOLDS = linspace(0, 9, 10)  # List of ints or floats.
THRESHOLDS = [2]
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

    omics = OmicsIntegration(rf"{DATA_PATH}/omics/raw_data_light.txt", samples_names={"SRR7984026Aligned.out.sam": "LL_1",
                                                                                      "SRR7984027Aligned.out.sam": "LL_2",
                                                                                      "SRR7984028Aligned.out.sam": "LL_3",
                                                                                      "SRR7984029Aligned.out.sam": "ML_1",
                                                                                      "SRR7984030Aligned.out.sam": "ML_2",
                                                                                      "SRR7984031Aligned.out.sam": "ML_3",
                                                                                      "SRR7984032Aligned.out.sam": "HL_1",
                                                                                      "SRR7984033Aligned.out.sam": "HL_2",
                                                                                      "SRR7984034Aligned.out.sam": "HL_3",
                                                                                      }, model=template_model)

    omics.getmm = read_csv(OMICS_DATA_PATH, index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics.sum_tech_reps()
    omics_data = omics.counts.applymap(lambda x: math.log2(x + 1))
    print(omics_data.describe())
    import matplotlib.pyplot as plt
    fig = omics_data.plot.density()
    fig.set_xticks([x / 2 for x in range(0, 40)])
    plt.savefig(rf"{DATA_PATH}/omics/counts_density_light.png")
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
                                                     thread_number=thred_number, omics_dataset=omics_data, params=params)

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

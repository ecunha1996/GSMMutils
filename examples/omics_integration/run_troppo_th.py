import math
import os
import pandas as pd
import sys

from GSMMutils.io import read_csv
from GSMMutils.omics.omics_integration import OmicsIntegration

sys.path.insert(0, r"/")
from GSMMutils.utils.config_variables import *
from GSMMutils.utils.pipeline_paths import *
from GSMMutils.omics.model_handle import load_model, sbml_model_reconstruction
from Tissue_specific_Reconstruction_Pipeline.Pipeline.troppo_integration import troppo_omics_integration
from Tissue_specific_Reconstruction_Pipeline.Pipeline.omics_processing import thresholding_filter
from Tissue_specific_Reconstruction_Pipeline.Pipeline.task_evaluation import task_eval


# TODO: Add the option for more integration algorithms in the pipeline.
# TODO: Add gap-filling to the pipeline.

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

    os.chdir(r"../")
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

    if THRESHOLDING_STRATEGY != 'default':
        omics_data = thresholding_filter(omics_dataframe=omics_data, thresholding_strategy=THRESHOLDING_STRATEGY,
                                         global_threshold_upper=GLOBAL_THRESHOLD_UPPER,
                                         global_threshold_lower=GLOBAL_THRESHOLD_LOWER,
                                         local_threshold=LOCAL_THRESHOLD)

        print(f'{THRESHOLDING_STRATEGY} threshold filter applied.')

    print('-------------------------------------------------------------------------------------------------------')
    print('------------------------------- Starting Omics Integration with Troppo. -------------------------------')
    print('-------------------------------------------------------------------------------------------------------')

    integration_result = {}

    for algorithm in ALGORITHMS:
        troppo_result_dict = {}
        print(f'Omics integration with {algorithm} started.')
        print('-------------------------------------------------------------------------------------------------------')

        if algorithm == 'fastcore':
            thred_number = THREAD_NUMBER_FASTCORE

        elif algorithm == 'tinit':
            thred_number = THREAD_NUMBER_TINIT

        else:
            return 'Algorithm not supported by this pipeline.'

        for threshold in INTEGRATION_THRESHOLDS:
            troppo_result = troppo_omics_integration(model=template_model, algorithm=algorithm, threshold=threshold,
                                                     thread_number=thred_number, omics_dataset=omics_data)

            for sample in list(troppo_result.keys()):
                th = str(round(threshold, 2)).replace('.', '_')
                troppo_result_dict[f'{MODEL}_{sample}_{algorithm}_t{th}'] = troppo_result[sample]
                integration_result[f'{MODEL}_{sample}_{algorithm}_t{th}'] = troppo_result[sample]

            print('----------------------------------------------------------'
                  '---------------------------------------------')

        if THRESHOLDING_STRATEGY == 'default':
            result_path = os.path.join(TROPPO_RESULTS_PATH,
                                       f'{MODEL}_{DATASET}_{algorithm}_{THRESHOLDING_STRATEGY}.csv')
        else:
            file_path = f'{MODEL}_{DATASET}_{algorithm}_{THRESHOLDING_STRATEGY}_' \
                        f'{GLOBAL_THRESHOLD_UPPER}_{GLOBAL_THRESHOLD_LOWER}_{LOCAL_THRESHOLD}.csv'
            result_path = os.path.join(TROPPO_RESULTS_PATH, file_path)

        troppo_result_dataframe = pd.DataFrame.from_dict(troppo_result_dict, orient='index')
        troppo_result_dataframe.to_csv(result_path)

    if EVALUATE_TASKS:
        print('------------------------ Starting task evaluation for the integration results. ------------------------')
        print('-------------------------------------------------------------------------------------------------------')

        task_eval(model_template=template_model, integration_results=integration_result)

        print('-------------------------------------------------------------------------------------------------------')

    if RECONSTRUCT_MODELS:
        print('----------------------- Starting reconstruction of the context-specific models. -----------------------')
        print('-------------------------------------------------------------------------------------------------------')

        for sample_name in list(integration_result.keys()):
            print(f'Context-specific model reconstruction for {sample_name} started.')
            print('---------------------------------------------------------------'
                  '----------------------------------------')

            sbml_model_reconstruction(model_template=template_model, sample=sample_name,
                                      integration_result_dict=integration_result)

            print('---------------------------------------------------------------'
                  '----------------------------------------')

    print('------------------------------------------ Pipeline Finished ------------------------------------------')
    print('-------------------------------------------------------------------------------------------------------')


if __name__ == '__main__':
    reconstruction_pipeline()

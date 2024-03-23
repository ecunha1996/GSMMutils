import traceback
import warnings
from os.path import abspath, dirname, join

from pandas import DataFrame

from ..utils.utils import block_print, enable_print
warnings.simplefilter("ignore")
CONFIG_PATH = abspath(join(dirname(__file__), '../../../config'))
block_print()
import cobra
import re
from troppo.omics.readers.generic import TabularReader
from troppo.methods_wrappers import ReconstructionWrapper
from cobamp.utilities.parallel import batch_run
enable_print()




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
        return [[k for k, v in data_map.get_scores().items() if (v is not None and v > threshold) or k in parameters['protected']]]

    def score_apply(data_map):
        dm = {k: 0 if v is None else (min(v, 10) - threshold) if k not in parameters['protected'] else 20
              for k, v in data_map.items()}
        return dm

    threshold, rec_wrapper, method = [parameters[parameter] for parameter in
                                      ['threshold', 'reconstruction_wrapper', 'algorithm']]

    # noinspection PyBroadException
    try:
        if method == 'fastcore':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=parameters['and_or_funcs'],
                                              integration_strategy=('custom', [integration_fx]), solver='CPLEX')

        elif method == 'tinit':
            return rec_wrapper.run_from_omics(omics_data=omics_container, algorithm=method, and_or_funcs=parameters['and_or_funcs'],
                                              integration_strategy=('continuous', score_apply), solver='CPLEX')

    except Exception as e:
        traceback.print_exc()
        return {r: False for r in rec_wrapper.model_reader.r_ids}


def troppo_omics_integration(model: cobra.Model, algorithm: str, threshold: float, thread_number: int,
                             omics_dataset: DataFrame, thresholds_map=None, params=None):
    """
    This function is used to run the Troppo's integration algorithms.

    Parameters
    ----------
    params
    thresholds_map
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
    protected = ["e_Biomass__cytop", ] + list(params['uptake_drains']['f2 medium'])
    parameters = {'threshold': threshold, 'reconstruction_wrapper': reconstruction_wrapper, 'algorithm': algorithm,
                  'protected': protected, 'and_or_funcs': (min, sum)
                  }

    batch_fastcore_res = batch_run(reconstruction_function, omics_data, parameters, threads=thread_number)

    result_dict = dict(zip([sample.condition for sample in omics_data], batch_fastcore_res))

    enable_print()
    print(f'Omics Integration with {details[1]} (Threshold = {details[2]}) Finished.')

    return result_dict

import os
from os.path import join

from .genome_annotation import GenomeAnnotation
from ..io import read_csv
from ..utils.remote import DockerClient


class InterProScan(GenomeAnnotation):
    def __init__(self):
        super().__init__()
        self._results = {}

    def load_results(self, path: str, name: str = None, **kwargs):
        """
        This function loads the results from the InterProScan.
        Parameters
        ----------
        name: str
            The name of the results.
        path: str
            The path to the InterProScan results.

        Returns
        -------

        """
        if not name:
            name = '.'.join(os.path.basename(path).split('.')[:-1])
        df = read_csv(path, sep="\t", header=None, **kwargs)
        df.columns = ['gene_id', 'hash', 'score', 'db', 'db_accession', 'db_description', 'start', 'end', 'evalue', 'status', 'date', 'interpro_accession', 'interpro_description']
        if name not in self.results:
            self.results[name] = {}
        self.results[name]['interproscan'] = df

    def load_results_from_folder(self, path: str):
        """
        This function loads the results from the InterProScan folder.
        Returns
        -------

        """
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(".tsv"):
                    self.load_results(join(root, file))


class InterProScanDocker(DockerClient):
    def __init__(self, data_directory=None, src_directory=None, examples_directory=None, config_directory=None, utilities_directory=None, server="turing",
                 interproscan_directory=None):
        super().__init__(data_directory, src_directory, examples_directory, config_directory, utilities_directory, server)
        self.interproscan_directory = interproscan_directory

    def build(self):
        cmd = f'{self.remote.container_tool} build {self.config_directory}/annotation -t interproscan'
        self.remote.run('interproscan', cmd)

    def run(self, interproscan_cmd):
        cmd = (
            f'{self.remote.container_tool} run --name interproscan --rm -v {self.data_directory}:/home/data:Z -v {self.src_directory}:/home/src/:Z '
            f'-v {self.examples_directory}:/home/examples/:Z -v {self.config_directory}:/home/config/:Z -v {self.utilities_directory}:/home/utilities/:Z -v {self.interproscan_directory}:/home/interproscan/:Z '
            f'interproscan sh -c "cd /home && /home/interproscan/interproscan.sh {interproscan_cmd}"'
        )
        self.remote.run('interproscan', cmd)

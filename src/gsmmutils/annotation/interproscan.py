import os
from os.path import join

from gsmmutils.utils.remote import Remote, DockerClient


class InterProScan:
    def __init__(self):
        pass


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

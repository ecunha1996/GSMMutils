import os
from os.path import join

from GSMMutils.utils.remote import Remote, DockerClient


class Busco:
    def __init__(self):
        pass


class BuscoDocker(DockerClient):
    def __init__(self, data_directory=None, src_directory=None, examples_directory=None, config_directory=None, utilities_directory=None, server="turing"):
        super().__init__(data_directory, src_directory, examples_directory, config_directory, utilities_directory, server)

    def pull(self):
        cmd = f'{self.remote.container_tool} pull ezlabgva/busco:v5.5.0_cv1'
        self.remote.run('busco', cmd)

    def run(self, busco_cmd):
        cmd = (
            f'{self.remote.container_tool} run -u $(id -u) --name busco --rm -v {self.data_directory}:/busco_wd:Z '
            f'ezlabgva/busco:v5.5.0_cv1 busco {busco_cmd}'
            )
        self.remote.run('busco', cmd)

import os
from os.path import join

from GSMMutils.utils.remote import Remote, DockerClient


class InterProScan:
    def __init__(self):
        pass


class InterProScanDocker(DockerClient):
    def __init__(self, data_directory=None, src_directory=None, examples_directory=None, config_directory=None, utilities_directory=None, server="turing"):
        super().__init__(data_directory, src_directory, examples_directory, config_directory, utilities_directory, server)

    def build(self):
        cmd = f'podman build {self.config_directory}/annotation -t interproscan'
        self.remote.run('interproscan', cmd)

    def run(self):
        cmd = (
            f'podman run --name interproscan -v {self.data_directory}:/home/data:Z -v {self.src_directory}:/home/src/:Z '
            f'-v {self.examples_directory}:/home/examples/:Z -v {self.config_directory}:/home/config/:Z '
            f'interproscan sh -c "python3 /home/examples/omics_integration/run_troppo_th.py && python3 /home/examples/sampling.py"'
            )
        self.remote.run('interproscan', cmd)

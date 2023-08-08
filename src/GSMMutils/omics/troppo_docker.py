import os
from os.path import join

from GSMMutils.utils.remote import Remote


class TroppoDocker:
    def __init__(self, data_directory=None, src_directory=None, examples_directory=None, config_directory=None):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("data")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.examples_directory = examples_directory or join(os.path.dirname(os.path.abspath(__file__)).split("examples")[0], "examples")
        self.config_directory = config_directory or join(os.path.dirname(os.path.abspath(__file__)).split("config")[0], "config")
        self.remote = Remote(data_directory, src_directory)

    def build(self):
        cmd = f'podman build . -t omics'
        self.remote.run('omics', cmd)

    def run(self):
        cmd = (f'podman run --name omics -v {self.data_directory}:/home/data:Z -v {self.src_directory}:/home/src/:Z -v {self.examples_directory}:/home/examples/:Z -v {self.config_directory}:/home/config/:Z '
               f'omics sh -c "python3 /home/examples/omics_integration/run_troppo_th.py && python3 /home/examples/sampling.py"'
               )  # > /home/data/output.txt
        self.remote.run('omics', cmd)

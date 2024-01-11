import os
from os.path import join
from gsmmutils.utils.remote import Remote


class OptimizationDocker:
    def __init__(self, data_directory=None, src_directory=None, examples_directory=None, config_directory=None, utilities_directory=None):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.examples_directory = examples_directory or join(
            os.path.dirname(os.path.abspath(__file__)).split("examples")[0], "examples")
        self.config_directory = config_directory or join(os.path.dirname(os.path.abspath(__file__)).split("config")[0],
                                                         "config")
        self.utilities_directory = utilities_directory or join(os.path.dirname(os.path.abspath(__file__)).split("utilities")[0], "utilities")
        self.remote = Remote(data_directory, src_directory)

    def build(self):
        cmd = f'podman build {self.config_directory}/optimization -t optimization'
        self.remote.run('optimization', cmd)

    def run(self):
        cmd = (f'podman run --name optimization --shm-size=10.24gb -v {self.data_directory}:/home/data/:Z -v '
               f'{self.src_directory}:/home/src/:Z -v {self.examples_directory}:/home/examples/:Z -v {self.config_directory}:/home/config/:Z '
               f'optimization '
               f'sh -c "python3 /home/examples/optimization/optimization.py > /home/data/output.txt"')  #
        self.remote.run('optimization', cmd)

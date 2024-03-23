import os
from os.path import join
from ..utils.remote import Remote


class DFBA:
    def __init__(self, data_directory=None, src_directory=None):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.remote = Remote(data_directory, src_directory)

    def run(self):
        cmd = (f'podman run --name dfba2 -v {self.data_directory}:/home/data/:Z -v '
               f'{self.src_directory}:/home/src/:Z davidtourigny/dfba sh -c "/usr/bin/python3 -m pip install --upgrade pip && pip install -q --root-user-action=ignore matplotlib seaborn scipy '
               f'joblib openpyxl psutil plotly kaleido timeout-decorator tqdm parallelbar && '
               f'python3 /home/src/gsmmutils/dynamic/run_dfba.py > /home/data/output.txt"')  #
        self.remote.run('dfba2', cmd)

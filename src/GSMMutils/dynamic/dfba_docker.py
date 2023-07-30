import os
from os.path import join
from GSMMutils.utils.remote import Remote


class DFBA:
    def __init__(self, data_directory=None, src_directory=None):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.remote = Remote(data_directory, src_directory)

    def run(self):
        cmd = (f'podman run --name dfba -v {self.data_directory}:/home/data -v '
               f'{self.src_directory}:/home/src/ davidtourigny/dfba sh -c "pip install -q matplotlib seaborn scipy '
               f'joblib openpyxl psutil plotly kaleido timeout-decorator tqdm parallelbar && python3 '
               f'/home/src/GSMMutils/dynamic/run_dfba.py > /home/data/output.txt"')  #
        self.remote.run('dfba', cmd)

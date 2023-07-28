import json
import os
from os.path import join
import paramiko

from GSMMutils.utils.remote import Remote
from GSMMutils.utils.utils import get_login_info


class dFBA:
    def __init__(self, data_directory=None, src_directory=None):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.remote = Remote(data_directory, src_directory)

    def run(self):
        cmd = f'podman run --name dfba -v {self.data_directory}:/home/data -v {self.src_directory}:/home/src/ davidtourigny/dfba sh -c "pip install -q matplotlib seaborn scipy joblib openpyxl psutil plotly kaleido timeout-decorator tqdm parallelbar && python3 /home/src/GSMMutils/dynamic/run_dfba.py > /home/data/output.txt"'  #
        self.remote.run('dfba', cmd)

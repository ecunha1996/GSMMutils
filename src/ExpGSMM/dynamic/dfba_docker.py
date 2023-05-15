import json
import os
from os.path import join
from pprint import pprint

import paramiko



class dFBA:
    def __init__(self, data_directory=None, src_directory=None):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.client = None
        self.login()


    def login(self):
        self.client = paramiko.SSHClient()
        self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.client.connect('turing.di.uminho.pt', username='ecunha', password='Somosporto1996')

    def run(self):
        stdin, stdout, stderr = self.client.exec_command('podman ps -a --format json')
        containers = json.loads(stdout.read().decode('utf-8'))
        exists = any([container['Names'][0] == 'dfba' for container in containers])
        if exists:
            print("Container already exists, running...")
            stdin, stdout, stderr = self.client.exec_command('podman start dfba')
        else:
            cmd = f'podman run --name dfba -v {self.data_directory}:/home/data -v {self.src_directory}:/home/src/ davidtourigny/dfba sh -c "pip install -q matplotlib seaborn scipy joblib openpyxl psutil plotly kaleido timeout-decorator tqdm parallelbar && python3 /home/src/ExpGSMM/dynamic/run_dfba.py > /home/data/output.txt"' #
            print(cmd)
            stdin, stdout, stderr = self.client.exec_command(cmd)

        print(stdout.read().decode('utf-8'))
        print(stderr.read().decode('utf-8'))

        # stdin, stdout, stderr = self.client.exec_command('podman rm dfba')

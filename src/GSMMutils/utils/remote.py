import json
import os
from os.path import join

import paramiko

from GSMMutils.utils.utils import get_login_info


class Remote:
    def __init__(self, data_directory=None, src_directory=None):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.client = None
        self.login()

    def login(self):
        try:
            self.client = paramiko.SSHClient()
            self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            host, username, password = get_login_info()
            self.client.connect(host, username=username, password=password)
            print(f"Connected to remote server {host}")
        except Exception as e:
            print(e)

    def run(self, docker_name, docker_cmd):
        _, stdout, stderr = self.client.exec_command('podman ps -a --format json')
        containers = json.loads(stdout.read().decode('utf-8'))
        exists = any([container['Names'][0] == docker_name for container in containers])
        if exists:
            print("Container already exists, running...")
            _, stdout, stderr = self.client.exec_command(f'podman start {docker_name}')
        else:
            print(docker_cmd)
            _, stdout, stderr = self.client.exec_command(docker_cmd)

        print(stdout.read().decode('utf-8'))
        print(stderr.read().decode('utf-8'))

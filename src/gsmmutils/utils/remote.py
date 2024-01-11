import json
import os
from os.path import join
import paramiko
from gsmmutils.utils.utils import get_login_info


class Remote:
    def __init__(self, data_directory=None, src_directory=None, server="turing"):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.client, self.container_tool  = None, None
        self.server = server
        self.login()

    def login(self):
        try:
            self.client = paramiko.SSHClient()
            self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            host, username, password, container_tool = get_login_info(self.server)
            self.client.connect(host, username=username, password=password)
            self.container_tool = container_tool
            print(f"Connected to remote server {host}")
        except Exception as e:
            print(e)

    def run(self, docker_name, docker_cmd):
        try:
            _, stdout, stderr = self.client.exec_command(f'{self.container_tool} ps -a --format json')
            stdout_as_str = stdout.read().decode('utf-8')
            containers = json.loads(stdout.read().decode('utf-8'))
            exists = any([container['Names'][0] == docker_name for container in containers])
        except:
            exists = False
        if exists:
            print("Container already exists, running...")
            _, stdout, stderr = self.client.exec_command(f'{self.container_tool} start {docker_name}')
        else:
            print(docker_cmd)
            _, stdout, stderr = self.client.exec_command(docker_cmd)

        print(stdout.read().decode('utf-8'))
        print(stderr.read().decode('utf-8'))


class DockerClient:
    def __init__(self, data_directory=None, src_directory=None, examples_directory=None, config_directory=None, utilities_directory=None, server="turing"):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("data")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.examples_directory = examples_directory or join(os.path.dirname(os.path.abspath(__file__)).split("examples")[0], "examples")
        self.config_directory = config_directory or join(os.path.dirname(os.path.abspath(__file__)).split("config")[0], "config")
        self.utilities_directory = utilities_directory or join(os.path.dirname(os.path.abspath(__file__)).split("utilities")[0], "utilities")
        self.server = server
        self.remote = Remote(data_directory, src_directory, self.server)



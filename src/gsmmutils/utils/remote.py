import json
import os
from os.path import join
import paramiko
from ..utils.utils import get_login_info
from stat import S_ISDIR, S_ISREG

class Remote:
    def __init__(self, data_directory=None, src_directory=None, server="turing"):
        """

        Parameters
        ----------
        data_directory
        src_directory
        server
        """
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.client, self.container_tool  = None, None
        self.server = server
        self.login()

    def login(self):
        """
        Method to login to remote server. It reads the login information from the config file and uses paramiko to connect to the server.
        If no login information is found, it will prompt the user to enter the login information.
        Returns
        -------

        """
        try:
            self.client = paramiko.SSHClient()
            self.client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            host, username, password, container_tool = get_login_info(self.server)
            self.client.connect(host, username=username, password=password)
            self.container_tool = container_tool
            print("Successfully connected to the remote server.")
        except Exception as e:
            print(e)

    def run(self, docker_name, docker_cmd):
        """
        Method to run a docker container on the remote server. If the container already exists, it will start the container.
        Parameters
        ----------
        docker_name
        docker_cmd

        Returns
        -------

        """
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

    def upload_data(self, localFilePath, remoteFilePath, isDirectory=False):
        """
        Method to upload data from local system to remote server.
        Parameters
        ----------
        localFilePath
        remoteFilePath

        Returns
        -------

        """
        sftp_client = self.client.open_sftp()
        try:
            sftp_client.put(localFilePath, remoteFilePath)
        except FileNotFoundError as err:
            print(f"File {localFilePath} was not found on the local system")
        sftp_client.close()

    def download_data(self, remoteFilePath, localFilePath, isDirectory=False):
        """
        Method to download data from remote server to local system.
        Parameters
        ----------
        remoteFilePath
        localFilePath

        Returns
        -------

        """
        sftp_client = self.client.open_sftp()
        try:
            if not isDirectory:
                sftp_client.get(remoteFilePath, localFilePath)
            else:
                sftp_get_recursive(remoteFilePath, localFilePath, sftp_client)
        except FileNotFoundError as err:
            print(f"File: {remoteFilePath} was not found on the source server {self.__hostName}:{self.__port}")
        finally:
            sftp_client.close()

    def sftp_put_recursive(self, path, dest, sftp):
        """
        Method to recursively upload data from local system to remote server.
        Parameters
        ----------
        path
        dest
        sftp

        Returns
        -------

        """
        if os.path.isdir(path):
            try:
                sftp.mkdir(dest)
            except IOError:
                pass
            for f in os.listdir(path):
                sftp_put_recursive(os.path.join(path, f), os.path.join(dest, f), sftp)
        else:
            sftp.put(path, dest)

    def sftp_get_recursive(self, path, dest, sftp):
        """
        Method to recursively download data from remote server to local system.
        Parameters
        ----------
        path
        dest
        sftp

        Returns
        -------

        """
        item_list = sftp.listdir_attr(path)
        dest = str(dest)
        if not os.path.isdir(dest):
            os.makedirs(dest, exist_ok=True)
        for item in item_list:
            mode = item.st_mode
            if S_ISDIR(mode):
                sftp_get_recursive(path + "/" + item.filename, dest + "/" + item.filename, sftp)
            else:
                sftp.get(path + "/" + item.filename, dest + "/" + item.filename)

    def exec(self, cmd):
        """
        Method to execute a command on the remote server.
        Parameters
        ----------
        cmd

        Returns
        -------

        """
        try:
            self.client.exec_command(cmd)
        except Exception as e:
            print(e)


class DockerClient:
    def __init__(self, data_directory=None, src_directory=None, examples_directory=None, config_directory=None, utilities_directory=None, server="turing"):
        self.data_directory = data_directory or join(os.path.dirname(os.path.abspath(__file__)).split("data")[0], "data")
        self.src_directory = src_directory or join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src")
        self.examples_directory = examples_directory or join(os.path.dirname(os.path.abspath(__file__)).split("examples")[0], "examples")
        self.config_directory = config_directory or join(os.path.dirname(os.path.abspath(__file__)).split("config")[0], "config")
        self.utilities_directory = utilities_directory or join(os.path.dirname(os.path.abspath(__file__)).split("utilities")[0], "utilities")
        self.server = server
        self.remote = Remote(data_directory, src_directory, self.server)



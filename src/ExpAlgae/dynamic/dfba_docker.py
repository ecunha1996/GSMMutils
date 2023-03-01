import os
from os.path import join

import docker
client = docker.from_env()



class dFBA:
    def __init__(self, data_directory):
        self.data_directory = os.getcwd()


    def run(self):
        containers = client.containers.list(all=True)
        exists = any([container.name == 'dfba' for container in containers])
        if exists:
            container = client.containers.get('dfba')
            container.exec_run('cd /home')
            exit_code, output = container.exec_run('python3 /home/src/ExpAlgae/dynamic/run_dfba.py', stdout=True, stderr=True)
        else:
            container = client.containers.run('davidtourigny/dfba',
                                              volumes={self.data_directory: {'bind': '/home/data/', 'mode': 'rw'}, join(os.path.dirname(os.path.abspath(__file__)).split("src")[0], "src"): {'bind': '/home/src/', 'mode': 'rw'}},
                                              tty=True,
                                              detach=True,
                                              name='dfba')
            container.exec_run('cd /home')
            container.exec_run('pip install scipy')
            container.exec_run('pip install matplotlib')
            container.exec_run('pip install openpyxl')
            container.exec_run('pip install -U kaleido')
            container.exec_run('pip install -U plotly')
            exit_code, output = container.exec_run('python3 /home/src/ExpAlgae/dynamic/run_dfba.py', stdout=True, stderr=True)
            container.stop()
        if exit_code: print(exit_code)
        if output: print(output)
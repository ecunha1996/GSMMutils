import time
from os.path import join

import matlab.engine

from ExpGSMM import DATA_PATH, DFBALAB_PATH


class DAE:
    def __init__(self, data_directory=join(DATA_PATH, "dynamic/DAE"), start:bool=True):
        self.data_directory = data_directory
        self.eng = None
        if start:
            self.start()

    def start(self):
        """
        Starts the matlab engine and adds the gurobi path to the matlab engine
        :return:
        """
        self.eng = matlab.engine.start_matlab()
        self.eng.cd(DFBALAB_PATH)
        self.eng.addpath(r'C:\gurobi912\win64\matlab')
        print("Matlab engine started")


    def run(self):
        """
        Runs the main function of the matlab engine
        :return:
        """
        self.eng.main("", 29)


    def batch_run(self):
        """
        Runs the main function of the matlab engine
        :return:
        """
        self.eng.AllTrials()
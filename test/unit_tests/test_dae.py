
import unittest
import os
import matlab.engine
from ExpAlgae.dynamic.DAE import DAE

class TestDAE(unittest.TestCase):

    def setUp(self):
        self.dae = DAE()

    def test_data_directory(self):
        self.assertEqual(type(self.dae.data_directory), str)

    def test_start_matlab_engine(self):
        self.dae.start()
        self.assertIsInstance(self.dae.eng, matlab.matlabengine.MatlabEngine)

    def test_run_main_function(self):
        self.dae.start()
        self.dae.run()


    def tearDown(self):
        self.dae = None

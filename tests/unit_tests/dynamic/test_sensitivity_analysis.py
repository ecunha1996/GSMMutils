import unittest
from os import chdir
from os.path import join, dirname
from gsmmutils.dynamic.sensitivity_analysis import SensitivityAnalysis
from gsmmutils.experimental.exp_matrix import ExpMatrix
from gsmmutils.graphics.plot import lineplot


class TestSensitivityAnalysis(unittest.TestCase):
    def setUp(self):
        data_path = join(dirname(__file__), '../../data')
        chdir(data_path)
        path = "Matrix_test.xlsx"
        self.matrix = ExpMatrix(path, conditions="Resume")
        self.analysis = SensitivityAnalysis(self.matrix)
        self.analysis.load_parameters(f"initial_parameters.json")

    def test_evaluate_dynamic_expression(self):
        # res = self.analysis.evaluate_dynamic_expression("starch_production", "Nitrogen_quota")
        # assert len(set(res.values())) == 1
        #
        # res = self.analysis.evaluate_dynamic_expression("starch_production", "x_storage")
        # assert len(set(res.values())) > 1
        #
        # res = self.analysis.evaluate_dynamic_expression("tag", "n_quota")
        # assert len(set(res.values())) > 1
        #
        # res = self.analysis.evaluate_dynamic_expression("chlorophyll", "p_quota")
        # lineplot(list(res.keys()), list(res.values()), xlabel="p_quota", ylabel="chlorophyll")
        #
        # res = self.analysis.evaluate_dynamic_expression("chlorophyll", "n_quota", param_range=[4, 7, 10])
        # graph = lineplot(list(res.keys()), list(res.values()), xlabel="n_quota", ylabel="chlorophyll")
        # graph.axhline()

        res = self.analysis.evaluate_dynamic_expression("carotene", "n_quota")
        lineplot(list(res.keys()), list(res.values()), xlabel="n_quota", ylabel="carotene")


if __name__ == '__main__':
    unittest.main()

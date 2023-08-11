import unittest

import pandas as pd
from statsmodels.multivariate.multivariate_ols import MultivariateTestResults

from GSMMutils.stats.stats import StatisticalAnalysis


class TestStats(unittest.TestCase):

    def setUp(self) -> None:
        self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/ancova/manova_data.csv")
        self.stats = StatisticalAnalysis(self.data)

    def test_anova(self):
        self.assertEqual(self.stats.anova('height ~ plant_var'), (None, None))

    def test_get_correlation(self):
        self.assertEqual(self.stats.get_correlation(), self.data.corr(method="pearson"))
        self.assertEqual(self.stats.get_correlation(columns=["height", "canopy_vol"]),
                         self.data[["height", "canopy_vol"]].corr(method="pearson"))
        self.assertEqual(self.stats.get_correlation(columns="height"), self.data["height"].corr(method="pearson"))
        self.assertEqual(self.stats.get_correlation(columns=["height", "canopy_vol"], method="kendall"),
                         self.data[["height", "canopy_vol"]].corr(method="kendall"))
        self.assertEqual(self.stats.get_correlation(columns="height", method="kendall"),
                         self.data["height"].corr(method="kendall"))
        self.assertEqual(self.stats.get_correlation(columns=["height", "canopy_vol"], method="spearman"),
                         self.data[["height", "canopy_vol"]].corr(method="spearman"))
        self.assertEqual(self.stats.get_correlation(columns="height", method="spearman"),
                         self.data["height"].corr(method="spearman"))

    def test_manova(self):
        self.assertIsInstance(self.stats.manova('height + canopy_vol ~ plant_var'), MultivariateTestResults)


if __name__ == "__main__":
    unittest.main()

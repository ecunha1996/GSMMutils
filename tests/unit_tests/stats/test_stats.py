import unittest

import pandas as pd
import pandas.testing as pd_testing
from statsmodels.multivariate.multivariate_ols import MultivariateTestResults
from statsmodels.regression.linear_model import RegressionResultsWrapper

from GSMMutils.stats.stats import StatisticalAnalysis


class TestStats(unittest.TestCase):
    def assertDataframeEqual(self, a, b):
        try:
            pd_testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException("DataFrames are not equal") from e

    def setUp(self) -> None:
        self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/ancova/manova_data.csv")
        self.stats = StatisticalAnalysis(self.data)

    def test_anova(self):
        anova_table, model = self.stats.anova('height ~ plant_var')
        self.assertEqual(anova_table.shape, (2, 4))
        self.assertIsInstance(model, RegressionResultsWrapper)

    def test_get_correlation(self):
        self.assertRaises(ValueError, self.stats.get_correlation)
        self.assertDataframeEqual(self.stats.get_correlation(cols=["height", "canopy_vol"]),
                                  self.data[["height", "canopy_vol"]].corr(method="pearson").round(2))
        self.assertDataframeEqual(self.stats.get_correlation(cols="height"), self.data[["height"]].corr(method="pearson"))
        self.assertDataframeEqual(self.stats.get_correlation(cols=["height", "canopy_vol"], method="kendall"),
                                  self.data[["height", "canopy_vol"]].corr(method="kendall").round(2)),
        self.assertDataframeEqual(self.stats.get_correlation(cols="height", method="kendall"),
                                  self.data[["height"]].corr(method="kendall").round(2)),
        self.assertDataframeEqual(self.stats.get_correlation(cols=["height", "canopy_vol"], method="spearman"),
                                  self.data[["height", "canopy_vol"]].corr(method="spearman").round(2)),
        self.assertDataframeEqual(self.stats.get_correlation(cols="height", method="spearman"),
                                  self.data[["height"]].corr(method="spearman").round(2))

    def test_manova(self):
        self.assertIsInstance(self.stats.manova('height + canopy_vol ~ plant_var'), MultivariateTestResults)


if __name__ == "__main__":
    unittest.main()

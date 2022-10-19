from statsmodels.formula.api import ols
from statsmodels.multivariate.manova import MANOVA
import statsmodels.api as sm


class StatisticalAnalysis:
    def __init__(self, data):
        self.data = data

    def get_correlation(self, columns: list = "all", method: str = "pearson"):
        if columns == "all":
            return self.data.corr(method=method)
        else:
            return self.data[columns].corr(method=method)

    def anova(self, formula):
        model = ols(formula, data=self.data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        return anova_table, model


    def manova(self, formula):
        fit= MANOVA.from_formula(formula, data=self.data)
        return fit.mv_test()


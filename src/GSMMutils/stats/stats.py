import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.multivariate.manova import MANOVA


def print_table(anova_table, title):
    title = title.center(45)
    print(f" {title} \n {anova_table} \n {'#' * 50}")


class StatisticalAnalysis:
    def __init__(self, data):
        self.data = data

    def get_correlation(self, columns: list | str = "all", method: str = "pearson"):
        if columns == "all":
            return self.data.corr(method=method)
        else:
            return self.data[columns].corr(method=method)

    def anova(self, formula, to_print=True):
        model = ols(formula, data=self.data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        if to_print:
            print_table(anova_table, formula)
        return anova_table, model

    def manova(self, formula):
        fit = MANOVA.from_formula(formula, data=self.data)
        test = fit.mv_test()
        print_table(test, formula)
        return test

import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.multivariate.manova import MANOVA


def print_table(anova_table, title):
    title = title.center(45)
    print(f" {title} \n {anova_table} \n {'#' * 50}")


class StatisticalAnalysis:
    def __init__(self, data):
        self.data = data

    def get_correlation(self, cols: list | str = "all", method: str = "pearson"):
        if cols == "all":
            if any(self.data.dtypes == "object"):
                raise ValueError("Cannot calculate correlation for non-numeric columns")
            return self.data.corr(method=method).round(2)
        else:
            if isinstance(cols, str):
                cols = [cols]
            if any(e == "object" for e in cols):
                raise ValueError("Cannot calculate correlation for non-numeric columns")
            return self.data[cols].corr(method=method).round(2)

    def anova(self, formula: str, to_print: bool = True) -> tuple:
        """
        Function to perform ANOVA test
        Parameters
        ----------
        formula: str
            formula for ANOVA
        to_print: bool
            print the table or not

        Returns
        -------
        anova_table: DataFrame
        model: OLS
        """
        model = ols(formula, data=self.data).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        if to_print:
            print_table(anova_table, formula)
        return anova_table, model

    def manova(self, formula):
        """
        Parameters
        ----------
        formula: str
            formula for MANOVA

        Returns
        -------
        test: MANOVA
        """
        fit = MANOVA.from_formula(formula, data=self.data)
        test = fit.mv_test()
        print_table(test, formula)
        return test

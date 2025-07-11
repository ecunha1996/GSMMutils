from pprint import pprint

import pandas as pd
from matplotlib import pyplot as plt

from gsmmutils import DATA_PATH

pd.set_option('display.max_columns', None)
from gsmmutils.experimental.biomass import Biomass
from gsmmutils.experimental.exp_matrix import ExpMatrix
from gsmmutils.graphics.plot import boxplot, qqplot
from gsmmutils.stats.stats import StatisticalAnalysis


def stats(matrix):
    matrix.conditions = matrix.conditions.rename({"[N] mmol": "N", "[P] mmol": "P", "Salinity g/L": "salinity", "Aeration rate": "aeration", 'growth_rate': 'umax', 'Productivity (g/L.h)': 'Pmax', 'Biomass (gDW/L)': 'biomass'}, axis=1)

    boxplot(matrix.conditions, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Pmax'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'Pmax': 'Pmax $(g \cdot L^{-1} \cdot d^{-1})$'})
    boxplot(matrix.conditions, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['biomass'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'biomass': 'Biomass $(g \cdot L^{-1} \cdot d^{-1})$'})
    boxplot(matrix.conditions, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['umax'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'umax': 'umax $(g \cdot L^{-1} \cdot d^{-1})$'})
    stats = StatisticalAnalysis(matrix.conditions)
    anova_table, model = stats.anova('biomass ~ P')
    print(anova_table)
    # hist(matrix.conditions, ['biomass'], title='Biomass', xlabel='Biomass $(g \cdot L^{-1})$', ylabel='Frequency')
    # hist(matrix.conditions, ['Pmax'], title='Maximum Productivity', xlabel='Pmax $(g \cdot L^{-1} \cdot h^{-1})$', ylabel='Frequency')
    # hist(matrix.conditions, ['umax'], title='Growth Rate', xlabel='Biomass $(h^{-1})$', ylabel='Frequency')
    qqplot(model, to_show=True)

if __name__ == '__main__':
    biomass = Biomass("e_Biomass__cytop", f"{DATA_PATH}/experimental/Biomass_exp_composition.xlsx")
    matrix = ExpMatrix(f"{DATA_PATH}/experimental/Matriz- DCCR Dunaliella salina_dfba_new.xlsx")
    matrix.conditions = "Resume"
    matrix.conditions = matrix.conditions.rename({"[N] mmol": "N", "[P] mmol": "P", "Salinity g/L": "salinity", "Aeration rate": "aeration", 'growth_rate': 'umax', 'Productivity (g/L.h)': 'Pmax', 'Biomass (gDW/L)': 'biomass'}, axis=1)
    stats(matrix)
    # m = pd.concat([biomass.biomass_matrix['macromolecules'], matrix.conditions[["N", "P", "salinity", "aeration"]]], axis=1)
    # m = m.dropna()
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Protein'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'Protein': 'Protein w/w'})
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Carbohydrate'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'Carbohydrate': 'Carbohydrate w/w'})
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Lipid'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'Lipid': 'Lipid w/w'})

    # results_dataframe = m
    # import scipy
    # slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(results_dataframe['salinity'], results_dataframe['Lipid'])
    #
    # # Create the regression line equation
    # equation = f'Lipid = {slope:.2f}salinity + {intercept:.2f}'
    #
    # # Create a scatter plot of the data points
    # plt.scatter(results_dataframe['salinity'], results_dataframe['Lipid'], label='Data Points')
    #
    # # Add the regression line to the plot
    # plt.plot(results_dataframe['salinity'], slope * results_dataframe['salinity'] + intercept, 'r', label='Regression Line')
    #
    # # Add labels and title to the plot
    # plt.xlabel('X-axis')
    # plt.ylabel('Y-axis')
    # plt.title('Linear Regression')
    # plt.legend()
    # print(equation)
    # # Display the equation on the plot
    # plt.text(1, 4, equation, fontsize=12)
    #
    # # Show the plot
    # plt.show()
    #
    #
    #
    #
    # stats = StatisticalAnalysis(m)
    # cor = stats.get_correlation()
    #
    # pprint(cor)
    #
    # m = pd.concat([biomass.biomass_matrix['pigments'].filter(regex='.*mean.*',axis=1), matrix.conditions[["N", "P", "salinity", "aeration"]]], axis=1)
    # m = m.dropna()
    # m.rename(columns = {"B-carotene (mean)": 'carotene', "Lutein (mean)": 'lutein', "Chlorophyll a (mean)": 'chla', "Chlorophyll b (mean)": 'chlb'}, inplace=True)
    #
    # m['caro_lutein'] = m['carotene']/m['lutein']
    #
    #
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['carotene'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'carotene': 'carotene w/w'})
    #
    #
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['chla'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'chla': 'chla w/w'})
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['chlb'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'chlb': 'chlb w/w'})
    #
    # m['lutein_concentration'] = m['lutein'] * matrix.conditions['biomass']
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['lutein'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'lutein': 'lutein w/w'})
    #
    # boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['lutein_concentration'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
    #         , y_labels={'lutein_concentration': 'lutein g/L'})
    #
    # stats = StatisticalAnalysis(m)
    # cor = stats.get_correlation()
    #
    # print(cor)
    #
    # stats.anova('P ~ chla')
    # stats.anova('N ~ chla')
    # stats.anova('salinity ~ chla')
    # stats.anova('aeration ~ chla')
    #
    # stats.anova('P ~ chlb')
    # stats.anova('N ~ chlb')
    # stats.anova('salinity ~ chlb')
    # stats.anova('aeration ~ chlb')
    #
    # stats.anova('P ~ carotene')
    # stats.anova('N ~ carotene')
    # stats.anova('salinity ~ carotene')
    # stats.anova('aeration ~ carotene')
    #
    # stats.anova('P ~ lutein')
    # stats.anova('N ~ lutein')
    # stats.anova('salinity ~ lutein')
    # stats.anova('aeration ~ lutein')
    #
    # stats.anova('P ~ caro_lutein')
    # stats.anova('N ~ caro_lutein')
    # stats.anova('salinity ~ caro_lutein')
    # stats.anova('aeration ~ caro_lutein')
    #
    # stats.manova('salinity + aeration ~ carotene')
    # stats.manova('P + salinity ~ carotene')
    # stats.manova('N + salinity ~ carotene')
    # stats.manova('N + P ~ carotene')

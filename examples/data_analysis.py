from pprint import pprint

import pandas as pd

from ExpGSMM import DATA_PATH

pd.set_option('display.max_columns', None)
from ExpGSMM.experimental.Biomass import Biomass
from ExpGSMM.experimental.ExpMatrix import ExpMatrix
from ExpGSMM.graphics.plot import boxplot, qqplot
from ExpGSMM.stats.stats import StatisticalAnalysis

def stats(matrix):
    matrix.conditions = "Resume"
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
    matrix = ExpMatrix(f"{DATA_PATH}/experimental/Matriz- DCCR Dunaliella salina_new.xlsx")
    matrix.conditions = "Resume"
    matrix.conditions = matrix.conditions.rename({"[N] mmol": "N", "[P] mmol": "P", "Salinity g/L": "salinity", "Aeration rate": "aeration", 'growth_rate': 'umax', 'Productivity (g/L.h)': 'Pmax', 'Biomass (gDW/L)': 'biomass'}, axis=1)
    stats(matrix)
    m = pd.concat([biomass.biomass_matrix['macromolecules'], matrix.conditions[["N", "P", "salinity", "aeration"]]], axis=1)
    m = m.dropna()
    boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Protein'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'Protein': 'Protein w/w'})
    boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Carbohydrate'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'Carbohydrate': 'Carbohydrate w/w'})
    boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['Lipid'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'Lipid': 'Lipid w/w'})

    stats = StatisticalAnalysis(m)
    cor = stats.get_correlation()

    pprint(cor)

    m = pd.concat([biomass.biomass_matrix['pigments'].filter(regex='.*mean.*',axis=1), matrix.conditions[["N", "P", "salinity", "aeration"]]], axis=1)
    m = m.dropna()
    m.rename(columns = {"B-carotene (mean)": 'carotene', "Lutein (mean)": 'lutein', "Chlorophyll a (mean)": 'chla', "Chlorophyll b (mean)": 'chlb'}, inplace=True)



    boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['carotene'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'carotene': 'carotene w/w'})


    boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['chla'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'chla': 'chla w/w'})
    boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['chlb'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'chlb': 'chlb w/w'})
    boxplot(m, x_cols=['P', 'N', 'salinity', 'aeration'], y_cols=['lutein'], to_show=True, x_labels={'P': 'P (mM)', 'N': 'N (mM)', 'salinity': 'NaCl $(g \cdot L^{-1})$', 'aeration': 'aeration rate'}
            , y_labels={'lutein': 'lutein w/w'})

    stats = StatisticalAnalysis(m)
    cor = stats.get_correlation()

    print(cor)

    stats.anova('P ~ chla')
    stats.anova('N ~ chla')
    stats.anova('salinity ~ chla')
    stats.anova('aeration ~ chla')

    stats.anova('P ~ chlb')
    stats.anova('N ~ chlb')
    stats.anova('salinity ~ chlb')
    stats.anova('aeration ~ chlb')

    stats.anova('P ~ carotene')
    stats.anova('N ~ carotene')
    stats.anova('salinity ~ carotene')
    stats.anova('aeration ~ carotene')

    stats.anova('P ~ lutein')
    stats.anova('N ~ lutein')
    stats.anova('salinity ~ lutein')
    stats.anova('aeration ~ lutein')

    stats.manova('salinity + aeration ~ carotene')
    stats.manova('P + salinity ~ carotene')
    stats.manova('N + salinity ~ carotene')
    stats.manova('N + P ~ carotene')

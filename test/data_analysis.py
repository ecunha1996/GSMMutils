from pprint import pprint

import pandas as pd
pd.set_option('display.max_columns', None)
from ExpAlgae.experimental.Biomass import Biomass
from ExpAlgae.experimental.ExpMatrix import ExpMatrix
from ExpAlgae.graphics.plot import boxplot
from ExpAlgae.stats.stats import StatisticalAnalysis

if __name__ == '__main__':
    data_directory = r"../data"
    biomass = Biomass("e_Biomass__cytop", "Biomass_exp_composition.xlsx")
    matrix = ExpMatrix("Matriz- DCCR Dunaliella salina_new.xlsx")
    matrix.conditions = "Resume"
    matrix.conditions = matrix.conditions.rename({"[N] mmol": "N", "[P] mmol": "P", "Salinity g/L": "salinity", "Aeration rate": "aeration", 'growth_rate': 'umax', 'Productivity (g/L.h)': 'Pmax', 'Biomass (gDW/L)': 'biomass'}, axis=1)
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

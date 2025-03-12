import cobra
from gsmmutils.model import MyModel
from cobra.flux_analysis import pfba
import warnings
warnings.filterwarnings('ignore')
import logging
logging.getLogger('cobra').setLevel(logging.CRITICAL)


def QY_AQ(non_diel_model, diel_model):
    fba_sol_non_diel = pfba(non_diel_model).fluxes
    fba_sol_diel_model = pfba(diel_model).fluxes

    return fba_sol_non_diel, fba_sol_diel_model


def test():
    original_model = MyModel('Populus_iPop7188.xml', "BiomassRxn")
    # original_model = cobra.io.read_sbml_model('Populus_iPop7188.xml')
    original_model.reactions.ATPM_c.bounds = (0, 1000)
    non_producing = original_model.test_reaction("BiomassRxn")

    print(non_producing.loc[non_producing.Flux == 0])

    # original_model.objective = "BiomassRxn"
    # print(original_model.optimize().status)
    # print(original_model.optimize().objective_value)


if __name__ == '__main__':
    # test()
    # original_model = cobra.io.read_sbml_model('Populus_iPop7188.xml')
    diel_populus_model = cobra.io.read_sbml_model("diel_populus_model.xml")
    #
    # original_model.objective = "EX_light"
    # original_model.objective_direction = "max"
    diel_populus_model.objective = "EX_light_Day"
    diel_populus_model.objective_direction = "max"
    #
    # original_model.reactions.get_by_id("BiomassRxn").lower_bound = 0.01
    # original_model.reactions.get_by_id("BiomassRxn").upper_bound = 0.01
    #
    diel_populus_model.reactions.get_by_id("Biomass_Total").lower_bound = 0.01
    diel_populus_model.reactions.get_by_id("Biomass_Total").upper_bound = 0.01
    #
    # fba_sol_non_diel = original_model.optimize()
    fba_sol_diel_model = diel_populus_model.optimize()

    # solver = 'cplex'
    # original_model.solver = solver
    # diel_populus_model.solver = solver
    # print(original_model.solver.status)
    print(diel_populus_model.solver.status)
    # print(fba_sol_non_diel['EX_light'])
    print(fba_sol_diel_model['EX_light_Day'])

    # fba_sol_non_diel, fba_sol_diel_model = QY_AQ(original_model, diel_populus_model)

    # data_quantum_assimilation = {
    #     'Quantum Yield': [fba_sol_non_diel["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"] / - fba_sol_non_diel["EX_light"],
    #                       fba_sol_diel_model["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_Day"] / - fba_sol_diel_model["EX_light_Day"]],
    #
    #     'Assimilation Quotient': [fba_sol_non_diel["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"] / fba_sol_non_diel["PSII_RXN"],
    #                               fba_sol_diel_model["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_Day"] / fba_sol_diel_model["PSII_RXN_Day"]]}
    #
    # tabel = pd.DataFrame(data_quantum_assimilation)
    #
    # tabel.index = ["Original Model", "Created Diel Model"]
    #
    # tabel.to_csv('QY&AQ_Populus.csv', sep=',')
    #
    # print(tabel)

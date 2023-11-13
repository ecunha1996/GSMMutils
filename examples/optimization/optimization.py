from cobra.io import read_sbml_model
from mewpy import get_simulator
from mewpy.optimization import BPCY, WYIELD, EA
from mewpy.problems import GOUProblem


def main():
    model = read_sbml_model(rf"/home/data/models/iJN678.xml")
    model.objective = "Ec_biomass_SynMixo"
    model.exchanges.EX_photon_e.bounds = (-65, -65)
    model.exchanges.EX_co2_e.bounds = (-2.7, -2.7)
    model.exchanges.EX_hco3_e.bounds = (-1.1, -1.1)
    print(model.summary())
    BIOMASS_ID = "Ec_biomass_SynMixo"
    PRODUCT_ID = "EX_mal_DASH_L_e"
    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method='lMOMA')
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
    simul = get_simulator(model)
    essential_genes = simul.essential_genes()
    problem = GOUProblem(simul, fevaluation=[
        evaluator_1, evaluator_2], candidate_max_size=5, non_target=essential_genes)
    ea = EA(problem, max_generations=100, visualizer=True)
    final_pop = ea.run()
    df = ea.dataframe()
    df.to_csv(rf"/home/data/optimization/synecho_gou_mal.tsv", index=False, sep="\t")


if __name__ == '__main__':
    main()

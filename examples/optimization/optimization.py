import json

import pandas as pd
from cobra.io import read_sbml_model
from mewpy import get_simulator
from mewpy.optimization import BPCY, WYIELD, EA
from mewpy.problems import GOUProblem, GKOProblem
import platform

if platform.system() == "Linux":
    data_path = "/home/data/"
else:
    data_path = r"C:\Users\Bisbii\PythonProjects\gsmmutils\data"


def main():
    model = read_sbml_model(rf"{data_path}/models/iJN678.xml")
    model.objective = "Ec_biomass_SynMixo"
    model.exchanges.EX_photon_e.bounds = (-65, -65)
    model.exchanges.EX_co2_e.bounds = (-2.7, -2.7)
    model.exchanges.EX_hco3_e.bounds = (-1.1, -1.1)
    print(model.summary())
    BIOMASS_ID = "Ec_biomass_SynMixo"
    PRODUCT_ID = "EX_mal_DASH_L_e"
    evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method='lMOMA')
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID, min_biomass_percent=0.3)
    evaluator_2.method = 'lMOMA'
    simul = get_simulator(model)
    essential_genes = simul.essential_genes()
    problem = GKOProblem(simul, fevaluation=[evaluator_1,
                                             evaluator_2], candidate_max_size=5, non_target=essential_genes)
    ea = EA(problem, max_generations=100, visualizer=True)
    final_pop = ea.run()
    df = ea.dataframe()
    df.to_csv(rf"{data_path}/optimization/synecho_gou_mal.tsv", index=False, sep="\t")
    solution = final_pop[1]
    sim = problem.simulator
    res = sim.simulate(constraints=solution.constraints, method='ROOM')
    print(solution)
    print(res)
    with open(rf"{data_path}/optimization/synecho_gou_mal.txt", "w") as f:
        f.write(str(res.objective_value))

def test_solution():
    res = pd.read_csv(rf"{data_path}/optimization/synecho_gou_mal.tsv", sep="\t")
    model = read_sbml_model(rf"{data_path}/models/iJN678.xml")
    model.objective = "Ec_biomass_SynMixo"
    model.exchanges.EX_photon_e.bounds = (-65, -65)
    model.exchanges.EX_co2_e.bounds = (-2.7, -2.7)
    model.exchanges.EX_hco3_e.bounds = (-1.1, -1.1)
    BIOMASS_ID = "Ec_biomass_SynMixo"
    PRODUCT_ID = "EX_mal_DASH_L_e"
    # evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method='lMOMA')
    evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID, min_biomass_percent=0.3)
    simul = get_simulator(model)
    essential_genes = simul.essential_genes()
    problem = GOUProblem(simul, fevaluation=[
        evaluator_2], candidate_max_size=5, non_target=essential_genes)
    solutions = res['Modification']
    for json_string in solutions:
        solution = json.loads(json_string.replace("'", "\""))
        constr = problem.solution_to_constraints(solution)
        print(simul.simulate(method='lMOMA', constraints=constr, objective="Ec_biomass_SynMixo"))


if __name__ == '__main__':
    main()
    # test_solution()

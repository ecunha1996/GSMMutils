import json
from os import makedirs
from os.path import exists

import numpy as np
import pandas as pd
from cobra import Metabolite
from cobra.flux_analysis import pfba, flux_variability_analysis as fva, moma, room
from joblib import Parallel, delayed
from tqdm import tqdm

from .COBRAmodel import MyModel
from ..graphics.plot import basic_scatter


class FSEOF:
    """
    Class to perform the FSEOF analysis.

    Parameters
    ----------
    model: MyModel
        The model to be used.
    targets: list[str]
        The target metabolites.
    workdir: str
        The directory to save the results.
    """
    def __init__(self, model: MyModel, targets: list[str], workdir: str = "results/fseof"):
        self.workdir = workdir
        if not exists(workdir):
            makedirs(workdir)
        self.adjusted_model = None
        self.model = model
        self.targets: list[Metabolite] = [model.metabolites.get_by_id(target) for target in targets]
        self.in_bio_reactions = {}
        for target in self.targets:
            if not exists(f"{workdir}/{target.id}"):
                makedirs(f"{workdir}/{target.id}")
            self.model.create_demand(target.id)
            r = [reaction for reaction in target.reactions if reaction.id.startswith("e_")]
            assert len(r) == 1
            self.in_bio_reactions[target.id] = r[0]

    def run(self):
        old_solution = pfba(self.model)
        assert old_solution.status == 'optimal'
        assert old_solution[self.model.biomass_reaction] > 0
        self.remove_precursor_in_macro_reaction()
        for target in self.targets:
            target_original_flux = round(old_solution[self.in_bio_reactions[target.id].id] * abs(self.in_bio_reactions[target.id].metabolites[target]), 5)  # R03824 0.006688 - 0.004070 - 0.001474 - 0.000344  e-Pigment 0.003990
            carbon_uptake = round(old_solution['EX_C00011__dra'], 5)
            demand_id = f"DM_{target.id}"
            print(f"{target.name} flux: ", target_original_flux)
            self.adjusted_model.reactions.get_by_id(demand_id).bounds = (target_original_flux, target_original_flux)
            self.adjusted_model.reactions.get_by_id("EX_C00011__dra").bounds = (carbon_uptake, carbon_uptake)
            new_solution = pfba(self.adjusted_model)
            fva_sol = fva(self.adjusted_model, fraction_of_optimum=1.0, processes=6)
            assert new_solution.status == 'optimal'
            assert new_solution[self.adjusted_model.biomass_reaction] > 0
            final_results_as_df, complete_results, growth_rate = self.fseof(demand_id)
            basic_scatter(growth_rate, path=f"{self.workdir}/{target}/{target.id}_growth.png", title=f"{target.name} growth", xlabel=f"{target.name} production", ylabel="Growth rate")
            final_results_as_df.to_csv(f"{self.workdir}/{target.id}/fseof_simulation.csv")
            increased, decreased, maintained, increased_valid, decreased_valid = self.evaluate_reactions(final_results_as_df, demand_id, target)
            self.adjusted_model.reactions.get_by_id(demand_id).bounds = (target_original_flux, 1000)
            increased_valid = {k: v for k, v in increased_valid.items() if (fva_sol.loc[k]['minimum'] > 0) or (fva_sol.loc[k]['maximum'] < 0)}
            ordered_obj_fun = self.reaction_over_expression(increased_valid, 1.2, new_solution, demand_id, target)
            self.get_pathways(ordered_obj_fun, self.adjusted_model, target)
            self.get_genes_with_impact(ordered_obj_fun, self.adjusted_model, target)

    def adjust_macromolecule_in_biomass(self):
        model = self.model.copy()
        reaction = model.bio_reaction
        total = round(sum([abs(reaction.metabolites[reactant]) for reactant in reaction.reactants if reactant.id != "C00001__cyto" and reactant.id != "C00002__cyto"]), 3)
        assert total == 1
        st_target = reaction.metabolites[self.targets]
        reaction.add_metabolites({model.metabolites.get_by_id("M8419__cyto"): -st_target})
        total = round(sum([abs(reaction.metabolites[reactant]) for reactant in reaction.reactants if reactant.id != "C00001__cyto" and reactant.id != "C00002__cyto"]), 3)
        assert total == 1 - abs(st_target)
        for reactant in reaction.reactants:
            if reactant.id != "C00001__cyto" and reactant.id != "C00002__cyto":
                st = reaction.metabolites[reactant]
                reaction.add_metabolites({reactant: -st})
                reaction.add_metabolites({reactant: st / total})
        total = round(sum([abs(reaction.metabolites[reactant]) for reactant in reaction.reactants if reactant.id != "C00001__cyto" and reactant.id != "C00002__cyto"]), 3)
        assert total == 1
        return model

    def remove_precursor_in_macro_reaction(self):
        model = self.model.copy()
        already_done = set()
        in_bio_reactions = {key: model.reactions.get_by_id(value.id) for key, value in self.in_bio_reactions.items()}
        for metabolite, reaction in in_bio_reactions.items():
            if reaction not in already_done:
                macro_id = [e for e in reaction.products if e.id.startswith("e_")][0]
                complete_g_gdw = abs(model.bio_reaction.metabolites[macro_id])
                new_reaction = reaction.copy()
                new_reaction.id = f'{"__".join(new_reaction.id.split("__")[:-1])}_target__{new_reaction.id.split("__")[-1]}'
                model.add_reactions([new_reaction])
                reaction.bounds = (0, 0)
                g_gdw_map = {}
                for met in new_reaction.metabolites:
                    g_gdw_map[met.id] = new_reaction.metabolites[met] * met.formula_weight * complete_g_gdw / 1000

                complete_g_gdw -= sum(abs(st) for met, st in g_gdw_map.items() if met in [met.id for met in self.targets])

                for target in self.targets:
                    model.set_stoichiometry(new_reaction.id, target.id, 0)

                for met in new_reaction.reactants:
                    new_st = g_gdw_map[met.id] / complete_g_gdw / met.formula_weight * 1000
                    model.set_stoichiometry(new_reaction.id, met.id, new_st)
                already_done.add(reaction)
                model.set_stoichiometry(model.bio_reaction.id, macro_id.id, -round(complete_g_gdw, 5))  # 0.0163
        total = sum(abs(st) for met, st in model.bio_reaction.metabolites.items() if -1 < model.bio_reaction.metabolites[met] < 0)
        # normalize to one
        for reactant in model.bio_reaction.reactants:
            if -1 < model.bio_reaction.metabolites[reactant] < 0:
                model.set_stoichiometry(model.bio_reaction.id, reactant.id, model.bio_reaction.metabolites[reactant] / total)
        total = round(sum(abs(st) for met, st in model.bio_reaction.metabolites.items() if -1 < model.bio_reaction.metabolites[met] < 0), 5)
        assert total == 1
        self.adjusted_model = model

    def fseof(self, objective):
        self.adjusted_model.reactions.get_by_id(objective).bounds = (0, 1000)
        self.adjusted_model.reactions.get_by_id("EX_C00011__dra").bounds = (-1000, 1000)
        with self.adjusted_model as model:
            model.objective = objective
            maximum_target_production = model.optimize().objective_value
        print(f"Maximum {objective}: ", maximum_target_production)
        span = np.linspace(0.05, 0.95, 10).tolist()
        complete_results, growth = {}, {}
        with self.adjusted_model as model:
            model.objective = self.adjusted_model.biomass_reaction
            for i, value in enumerate(span):
                target_production = round(value * maximum_target_production, 4)
                model.reactions.get_by_id(objective).bounds = (target_production, target_production)
                try:
                    solution = pfba(model)
                    print(f"Percentage: {value:.2f}\t{objective}: {target_production}\tBiomass: {round(solution[model.biomass_reaction], 3)}")
                    complete_results[str(round(value, 3))] = solution
                    growth[target_production] = round(solution[model.biomass_reaction], 3)
                except Exception as e:
                    print(e)
        final_results_as_dict = {}
        for key, value in complete_results.items():
            final_results_as_dict[key] = value.fluxes
        final_results_as_df = pd.DataFrame.from_dict(final_results_as_dict, orient='index')
        return final_results_as_df, complete_results, growth

    def evaluate_reactions(self, final_results_as_df, objective, target):
        """
        Method to evaluate the reactions that are increased, decreased or maintained.
        Parameters
        ----------
        final_results_as_df: pd.DataFrame
            A DataFrame with the results of the FSEOF simulation.
        objective: str
            The id of the demand reaction.
        target: Metabolite
            The target metabolite.

        Returns
        -------
        increased: dict
            Dictionary with the reactions that are increased.
        decreased: dict
            Dictionary with the reactions that are decreased.
        maintained: dict
            Dictionary with the reactions that are maintained.
        increased_valid: dict
            Dictionary with the reactions that are increased and valid.
        decreased_valid: dict
            Dictionary with the reactions that are decreased and valid.
        """
        print("Evaluating reactions")
        increased, decreased, maintained, increased_valid, decreased_valid = {}, {}, {}, {}, {}
        reactions = final_results_as_df.columns
        for reaction in reactions:
            y = final_results_as_df[reaction].abs()
            x = final_results_as_df.index.astype(float)
            slope, intercept = np.polyfit(x, y, 1)
            slope = round(slope, 3)
            if slope > 0:
                increased[reaction] = slope
            elif slope < 0:
                decreased[reaction] = slope
            else:
                maintained[reaction] = slope
        span = np.linspace(0.05, 0.95, 10).tolist()
        fva_solution = {}
        with self.adjusted_model as temp_model:
            temp_model.objective = objective
            maximum_target_production = temp_model.optimize().objective_value
        for i in tqdm(range(len(span))):
            target_production = round(span[i] * maximum_target_production, 4)
            self.adjusted_model.reactions.get_by_id(objective).bounds = (target_production, target_production)
            fva_sol = fva(self.adjusted_model, list(increased.keys()) + list(decreased.keys()), fraction_of_optimum=1.0, processes=6)
            fva_solution[round(span[i], 3)] = fva_sol

        for reaction in increased.keys():
            valid_fva = self.validate_with_fva(fva_solution, reaction)
            if valid_fva:
                increased_valid[reaction] = increased[reaction]

        for reaction in decreased.keys():
            valid_fva = self.validate_with_fva(fva_solution, reaction)
            if valid_fva:
                decreased_valid[reaction] = decreased[reaction]

        print("Increased: ", len(increased))
        print("Decreased: ", len(decreased))
        print("Maintained: ", len(maintained))
        print("Increased valid: ", len(increased_valid))
        print("Decreased valid: ", len(decreased_valid))
        with open(f'{self.workdir}/{target.id}/increased.json', 'w') as fp:
            json.dump({k: v for k, v in sorted(increased.items(), key=lambda item: item[1])}, fp, indent=None, separators=(",\n", ":"))
        with open(f'{self.workdir}/{target.id}/decreased.json', 'w') as fp:
            json.dump({k: v for k, v in sorted(decreased.items(), key=lambda item: item[1])}, fp, indent=None, separators=(",\n", ":"))
        with open(f'{self.workdir}/{target.id}/maintained.json', 'w') as fp:
            json.dump({k: v for k, v in sorted(maintained.items(), key=lambda item: item[1])}, fp, indent=None, separators=(",\n", ":"))
        with open(f'{self.workdir}/{target.id}/increased_valid.json', 'w') as fp:
            json.dump({k: v for k, v in sorted(increased_valid.items(), key=lambda item: item[1])}, fp, indent=None, separators=(",\n", ":"))
        with open(f'{self.workdir}/{target.id}/decreased_valid.json', 'w') as fp:
            json.dump({k: v for k, v in sorted(decreased_valid.items(), key=lambda item: item[1])}, fp, indent=None, separators=(",\n", ":"))
        self.adjusted_model.get_reactions_pathways_map()
        pathways_counter = {}
        for reaction, value in increased_valid.items():
            pathways = self.adjusted_model.reactions_pathways_map[reaction]
            for pathway in pathways:
                if pathway in pathways_counter:
                    pathways_counter[pathway] += 1
                else:
                    pathways_counter[pathway] = 1

        with open(f'{self.workdir}/{target.id}/pathways_increased.json', 'w') as fp:
            json.dump(pathways_counter, fp, indent=4)
        pathways_counter = {}
        for reaction, value in decreased_valid.items():
            pathways = self.adjusted_model.reactions_pathways_map[reaction]
            for pathway in pathways:
                if pathway in pathways_counter:
                    pathways_counter[pathway] += 1
                else:
                    pathways_counter[pathway] = 1
        with open(f'{self.workdir}/{target.id}/pathways_decreased.json', 'w') as fp:
            json.dump(pathways_counter, fp, indent=4)
        return increased, decreased, maintained, increased_valid, decreased_valid

    @staticmethod
    def validate_with_fva(fva_solution, reaction):
        """

        Parameters
        ----------
        fva_solution
        reaction

        Returns
        -------

        """
        for key, solution in fva_solution.items():
            values = solution.loc[reaction]
            if (values["maximum"] > 0 and values["minimum"] > 0 and values["minimum"] > values["maximum"] / 2) or (values["maximum"] < 0 and values["minimum"] < 0 and values["minimum"] < values["maximum"] / 2):
                return True
        return False

    def reaction_over_expression(self, reaction_list, expression, original_simulation, demand_id, target) -> dict:
        """

        Parameters
        ----------
        reaction_list
        expression
        original_simulation
        demand_id

        Returns
        -------

        """
        print("Over expression")
        obj_fun = Parallel(n_jobs=5)(delayed(self.over_expression_in_parallel)(original_simulation, expression, demand_id, reaction) for reaction in reaction_list)
        obj_fun = {k: v for d in obj_fun for k, v in d.items()}
        ordered_obj_fun = {k: v for k, v in sorted(obj_fun.items(), key=lambda item: item[1], reverse=True)}
        print(ordered_obj_fun)
        with open(f'{self.workdir}/{target.id}/obj_fun.json', 'w') as fp:
            json.dump(ordered_obj_fun, fp, indent=4)
        return ordered_obj_fun

    def over_expression_in_parallel(self, original_simulation, expression, demand_id, reaction):
        obj_fun = {}
        with self.adjusted_model as model:
            try:
                model.objective = model.biomass_reaction
                original_flux = original_simulation[reaction]
                if round(original_flux, 5) != 0:
                    model.reactions.get_by_id(reaction).bounds = (round(expression * original_flux, 5), round(expression * original_flux, 5))
                    solution = moma(model, original_simulation, False)
                    if solution.status == 'optimal':
                        original_target = original_simulation[demand_id]
                        original_biomass = original_simulation[model.biomass_reaction]
                        moma_target_production = solution[demand_id]
                        moma_biomass = solution[model.biomass_reaction]
                        obj_fun[reaction] = (moma_target_production / original_target) * (moma_biomass / original_biomass)
                    else:
                        print(reaction)
                        print(solution.status)
                        print(model.reactions.get_by_id(reaction).bounds)
                        print("#" * 100)
                else:
                    print(reaction, " has flux equal to 0")
            except Exception as e:
                print(e)
                print("Reaction: ", reaction, " not feasible")
        return obj_fun

    def get_pathways(self, ordered_obj_fun, adjusted_model, target):
        """

        Parameters
        ----------
        ordered_obj_fun
        adjusted_model

        Returns
        -------

        """
        adjusted_model.get_reactions_pathways_map()
        pathways_counter = {}
        for reaction, value in ordered_obj_fun.items():
            if value > 0:
                pathways = adjusted_model.reactions_pathways_map[reaction]
                for pathway in pathways:
                    if pathway in pathways_counter:
                        pathways_counter[pathway] += 1
                    else:
                        pathways_counter[pathway] = 1
        print(pathways_counter)
        with open(f'{self.workdir}/{target.id}/pathways.json', 'w') as fp:
            json.dump(pathways_counter, fp, indent=4)

    def get_genes_with_impact(self, ordered_obj_fun, adjusted_model, target):
        genes_of_interest = {}
        for reaction, value in ordered_obj_fun.items():
            if value > 0:
                genes = adjusted_model.reactions.get_by_id(reaction).genes
                for gene in genes:
                    if 'ec-code' in adjusted_model.reactions.get_by_id(reaction).annotation:
                        genes_of_interest[gene.id] = adjusted_model.reactions.get_by_id(reaction).annotation['ec-code']
                    else:
                        genes_of_interest[gene.id] = "No EC code"
        print(genes_of_interest)
        with open(f'{self.workdir}/{target.id}/genes.txt', 'w') as fp:
            for gene in genes_of_interest.keys():
                fp.write(gene + "\n")
        with open(f'{self.workdir}/{target.id}/genes.json', 'w') as fp:
            json.dump(genes_of_interest, fp, indent=4)

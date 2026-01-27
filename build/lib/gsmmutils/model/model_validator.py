from os.path import join
from cobra.flux_analysis import pfba, find_blocked_reactions, flux_variability_analysis as fva, fastcc
from pandas import DataFrame
from .COBRAmodel import MyModel
from ..utils.configs import get_config
DATA_PATH = get_config().get("PATHS", "DATA_PATH")


import logging
logging.getLogger('cobra').setLevel(logging.CRITICAL)


class ModelValidator:
    def __init__(self, model: MyModel, atpm_reaction=None, cytoplasm_abb="__cytop"):
        self.old_model = None
        self.model = model
        self.atpm_reaction = atpm_reaction
        self.cytoplasm_abb = cytoplasm_abb

    def validate(self):
        self.check_reactions_equal_metabolites()
        self.check_balance()
        self.old_model = self.model.copy()
        self.check_energy_producing_cycles()
        self.model = self.old_model
        self.check_biomass_production()
        self.check_blocked_reactions()
        self.check_unbounded_flux()
        self.check_sbc()

    def check_balance(self, show_biomass_reactions=False):
        mass_balance = {}
        for reaction in self.model.reactions:
            balance = reaction.check_mass_balance()
            if (reaction.check_mass_balance() and "EX_" not in reaction.id and any([round(value, 5) != 0 for value in balance.values()])
                    and "DM_" not in reaction.id and "SK_" not in reaction.id):
                if show_biomass_reactions and str(reaction.id).startswith("e_"):
                    mass_balance[str(reaction.id)] = reaction.check_mass_balance()
                else:
                    if not str(reaction.id).startswith("e_"):
                        mass_balance[str(reaction.id)] = reaction.check_mass_balance()
        if mass_balance:
            print("Mass balance: NOT OK")
            for key in mass_balance:
                print(key, mass_balance[key])
        else:
            print("Mass balance: OK")
        return mass_balance

    def check_energy_producing_cycles(self):
        def check_energy_production_for_met(current_model: MyModel, reactant_product_pair: tuple[str, str],
                                            metabolites: dict):
            current_model.create_reaction(f"test_reaction_for_{reactant_product_pair[0]}",
                                          metabolites)
            for ex in current_model.exchanges:
                if ex.reactants[0].id not in to_keep:
                    ex.bounds = (0, 1000)
            current_model.objective = current_model.reactions.get_by_id(f"test_reaction_for_{reactant_product_pair[0]}")
            return current_model.optimize().objective_value

        old_bounds = (0, 1000)
        if self.atpm_reaction and self.atpm_reaction in [r.id for r in self.model.reactions]:
            old_bounds = self.model.reactions.get_by_id(self.atpm_reaction).bounds
            self.model.reactions.get_by_id(self.atpm_reaction).bounds = (0, 0)
        energy_cycles = {}
        energy_metabolites = [("C00002", "C00008"), ("C00044", "C00035"), ("C00063", "C00112"), ("C00075", "C00015"),
                              ("C00081", "C00104")]
        redox_metabolites = [("C01352", "C00016"), ("C00004", "C00003"), ("C00005", "C00006")]
        to_keep = ["C00001", "C00080", "C00009"]
        metabolites_in_model = {m.id for m in self.model.metabolites}
        for met in energy_metabolites:
            if {met[0]}.issubset(metabolites_in_model) and {met[0]}.issubset(metabolites_in_model):
                copy_model = self.model.copy()
                energy_cycles[met[0]] = check_energy_production_for_met(copy_model, met,
                                                                        {copy_model.metabolites.get_by_id(
                                                                            met[0] + self.cytoplasm_abb): -1,
                                                                         copy_model.metabolites.C00001__cytop: -1,
                                                                         copy_model.metabolites.get_by_id(
                                                                             met[1] + self.cytoplasm_abb): 1,
                                                                         copy_model.metabolites.C00009__cytop: 1,
                                                                         copy_model.metabolites.C00080__cytop: 1})

        for met in redox_metabolites:
            if {met[0]}.issubset(metabolites_in_model) and {met[0]}.issubset(metabolites_in_model):
                copy_model = self.model.copy()
                energy_cycles[met[0]] = check_energy_production_for_met(copy_model, met,
                                                                        {copy_model.metabolites.get_by_id(
                                                                            met[0] + self.cytoplasm_abb): -1,
                                                                         copy_model.metabolites.get_by_id(
                                                                             met[1] + self.cytoplasm_abb): 1,
                                                                         copy_model.metabolites.C00080__cytop: 1})
        if self.atpm_reaction and self.atpm_reaction in [r.id for r in self.model.reactions]:
            self.model.reactions.get_by_id(self.atpm_reaction).bounds = old_bounds
        if energy_cycles and any([energy_cycles[key] > 0 for key in energy_cycles]):
            print("Energy cycles: NOT OK")
            for key in energy_cycles:
                print(key, energy_cycles[key])
        else:
            print(f"Energy cycles: OK")
        return energy_cycles

    def check_biomass_production(self):
        sol = self.model.optimize()
        if sol.objective_value > 0:
            if sol.objective_value < 2:
                print(f"Biomass production OK: {sol.objective_value}")
            else:
                print(f"Biomass production might be too high: {sol.objective_value}")
        else:
            print("Biomass production is not satisfied!")

        pfba_sol = pfba(self.model)
        print(self.model.summary(pfba_sol))
        return pfba_sol

    def check_blocked_reactions(self):
        blocked = find_blocked_reactions(self.model, open_exchanges=True)
        with open(join(DATA_PATH, "blocked_reactions.txt"), "w") as f:
            for reaction in blocked:
                f.write(str(reaction) + "\n")
        print(f"Blocked reactions: {len(blocked)} -> ({round(len(blocked) / len(self.model.reactions), 2) * 100}%)")
        print("Running fastcc...")
        consistent_model = fastcc(self.model)
        print("Consistent reactions: ", len(consistent_model.reactions))
        print("Inconsistent reactions: ", len(self.model.reactions) - len(consistent_model.reactions))
        consistent_model_reaction_ids = set([r.id for r in consistent_model.reactions])
        return blocked, list(set(self.model.reaction_ids) - consistent_model_reaction_ids)

    def check_unbounded_flux(self) -> DataFrame:
        fva_sol = fva(self.model, fraction_of_optimum=1.0)
        sbc = fva_sol.loc[(fva_sol["minimum"].round(5) < 0) & (fva_sol["maximum"].round(5) > 0)]
        print(f"Reactions with Unbounded flux: {len(sbc)} -> ({round(len(sbc) / len(self.model.reactions), 2) * 100}%)")
        print(sbc)
        return sbc

    def check_sbc(self):
        for ex in self.model.exchanges:
            ex.bounds = (0, 0)
        if self.atpm_reaction and self.atpm_reaction in self.model.reaction_ids:
            self.model.reactions.get_by_id(self.atpm_reaction).bounds = (0, 1000)
        fva_sol = fva(self.model, fraction_of_optimum=1.0)
        sbc = fva_sol.loc[(fva_sol["minimum"].round(5) < 0) & (fva_sol["maximum"].round(5) > 0)]
        sbc.to_csv(join(DATA_PATH, "sbc.csv"))
        print(f"SBC: {len(sbc)} -> ({round(len(sbc) / len(self.model.reactions), 2) * 100}%)")
        return sbc

    def check_reactions_equal_metabolites(self):
        inconsistency = False
        exceptions = ["PRISM", "e_Pigment", "e_Biomass"]
        already_printed = []
        for reaction in self.model.reactions:
            for second_reaction in self.model.reactions:
                exception = False
                for e in exceptions:
                    if e in reaction.id or e in second_reaction.id:
                        exception = True
                if not exception:
                    if reaction.id != second_reaction.id:
                        if (reaction.reactants == second_reaction.reactants
                                and reaction.products == second_reaction.products):
                            string_message = f"Reactions {reaction.id} and {second_reaction.id} are equal!"
                            sorted_string_message = " ".join(sorted(string_message.split()))
                            if sorted_string_message not in already_printed:
                                already_printed.append(sorted_string_message)
                                print(string_message)
                            inconsistency = True
                        if (reaction.reactants == second_reaction.products and
                                reaction.products == second_reaction.reactants):
                            string_message = f"Reactions {reaction.id} and {second_reaction.id} are equal!"
                            sorted_string_message = " ".join(sorted(string_message.split()))
                            if sorted_string_message not in already_printed:
                                already_printed.append(sorted_string_message)
                                print(string_message)
                            inconsistency = True
        if not inconsistency:
            print("No duplicate reactions found!")
        return already_printed

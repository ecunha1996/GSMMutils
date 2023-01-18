
from cobra import Reaction
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba, find_blocked_reactions, flux_variability_analysis as fva


def load_model(model_path):
    print("Loading model...")
    model = read_sbml_model(model_path)
    print("Model name: {}".format(model.name) + "\n" + "Reactions: " + str(len(model.reactions)) + "\n" + "Metabolites: " + str(len(model.metabolites)) + "\n" + "Genes: " + str(len(model.genes)))
    return model


def define_media(model, media=None):
    for key in media.keys():
        model.exchanges.get_by_id(key).lower_bound = media[key]
    return model

def check_under_limit(reaction):
    balance= reaction.check_mass_balance()
    if balance:
        for key in balance:
            if round(balance[key],5) != 0:
                return False
    return True


def check_balance(model, show_biomass_reactions=False):
    mass_balance = {}
    for reaction in model.reactions:
        if reaction.check_mass_balance() and "EX_" not in reaction.id and not check_under_limit(reaction) and "DM_" not in reaction.id:
            if str(reaction.id.split('__')[0]) + str(reaction.check_mass_balance()) not in mass_balance:
                if str(reaction.id).startswith("e_"):
                    if show_biomass_reactions:
                        mass_balance[str(reaction.id)] =  reaction.check_mass_balance()
                else:
                    mass_balance[str(reaction.id)] = reaction.check_mass_balance()
    if mass_balance:
        print("Mass balance: NOT OK")
        for key in mass_balance:
            print(key, mass_balance[key])
    else:
        print("Mass balance: OK")
    return mass_balance


def check_energy_producing_cycles(model, atpm_reaction = "ATPm__cytop", cytoplasm_abb = "__cytop"):
    old_bounds = model.reactions.get_by_id(atpm_reaction).bounds
    model.reactions.get_by_id(atpm_reaction).bounds = (0, 0)
    energy_cycles = {}
    energy_metabolites = [("C00002", "C00008"),("C00044","C00035"),("C00063","C00112"),("C00075","C00015"),("C00081","C00104")]
    redox_metabolites = [("C01352","C00016"),("C00004","C00003") ,("C00005","C00006")]
    to_keep = ["C00001", "C00080", "C00009"]
    for met in energy_metabolites:
        copy_model = model.copy()
        new_reaction = Reaction(id=f"test_reaction_for_{met[0]}")
        new_reaction.add_metabolites({copy_model.metabolites.get_by_id(met[0]+ cytoplasm_abb):-1, copy_model.metabolites.C00001__cytop:-1, copy_model.metabolites.get_by_id(met[1]+cytoplasm_abb):1, copy_model.metabolites.C00009__cytop: 1, copy_model.metabolites.C00080__cytop: 1})
        copy_model.add_reaction(new_reaction)
        for ex in copy_model.exchanges:
            if ex.reactants[0].id not in to_keep:
                ex.bounds = (0, 1000)
        copy_model.objective = copy_model.reactions.get_by_id(f"test_reaction_for_{met[0]}")
        sol = copy_model.optimize()
        energy_cycles[met[0]] = sol.objective_value
    for met in redox_metabolites:
        copy_model = model.copy()
        new_reaction = Reaction(id=f"test_reaction_for_{met[0]}")
        new_reaction.add_metabolites({copy_model.metabolites.get_by_id(met[0]+ cytoplasm_abb):-1, copy_model.metabolites.get_by_id(met[1]+cytoplasm_abb):1, copy_model.metabolites.C00080__cytop: 1})
        copy_model.add_reaction(new_reaction)
        for ex in copy_model.exchanges:
            if ex.reactants[0].id not in to_keep:
                ex.bounds = (0, 1000)
        copy_model.objective = copy_model.reactions.get_by_id(f"test_reaction_for_{met[0]}")
        sol = copy_model.optimize()
        energy_cycles[met[0]] = sol.objective_value
    model.reactions.get_by_id(atpm_reaction).bounds = old_bounds
    if energy_cycles and any([energy_cycles[key] > 0 for key in energy_cycles]):
        print("Energy cycles: NOT OK")
        for key in energy_cycles:
            print(key, energy_cycles[key])
    else:
        print(f"Energy cycles: OK")
    return energy_cycles


def check_biomass_production(model):
    sol = model.optimize()
    if sol.objective_value > 0:
        if sol.objective_value < 2:
            print(f"Biomass production OK: {sol.objective_value}")
        else:
            print(f"Biomass production might be too high: {sol.objective_value}")
    else:
        print("Biomass production is not satisfied!")

    pfba_sol = pfba(model)
    print(model.summary(pfba_sol))


def check_blocked_reactions(model):
    blocked = find_blocked_reactions(model, open_exchanges=True)
    print(f"Blocked reactions: {len(blocked)} -> ({round(len(blocked)/len(model.reactions),2)*100}%)")


def check_unbounded_flux(model):
    fva_sol = fva(model, fraction_of_optimum=1.0)
    sbc = fva_sol.loc[(round(fva_sol["minimum"], 5) < 0) & (round(fva_sol["maximum"], 5) > 0)]
    print(f"Reactions with Unbounded flux: {len(sbc)} -> ({round(len(sbc)/len(model.reactions),2)*100}%)")
    print(sbc)


def check_sbc(model, atpm_reaction):
    for ex in model.exchanges: ex.bounds = (0, 0)
    model.reactions.get_by_id(atpm_reaction).bounds = (0,1000)
    fva_sol = fva(model, fraction_of_optimum=1.0)
    sbc = fva_sol.loc[(round(fva_sol["minimum"], 5) < 0) | (round(fva_sol["maximum"], 5) > 0)]
    print(f"SBC: {len(sbc)} -> ({round(len(sbc) / len(model.reactions), 2) * 100}%)")
    print(sbc)


def check_consistency(model, atpm_reaction, cytoplasm_abb = "__cytop"):
    check_balance(model)
    check_energy_producing_cycles(model, atpm_reaction = atpm_reaction, cytoplasm_abb = cytoplasm_abb)
    check_biomass_production(model)
    check_blocked_reactions(model)
    check_unbounded_flux(model)
    check_sbc(model, atpm_reaction = atpm_reaction)








from gsmmutils.model.COBRAmodel import *
from gsmmutils import DATA_PATH


def load():
    model = MyModel(join(DATA_PATH, 'models/model_pl.xml'), 'e_Biomass__in')
    return model


def filter_components_by_biomass(e_lipid_reactions):
    """

    Parameters
    ----------
    model
    e_lipid_reactions: dict
        Dictionary where the keys are the general lipids, and the value is a dict with the number of fatty acids of the lipid, and the reactants of the lipid.

    Returns
    -------

    """
    acyl_chains = set()
    for key, value in e_lipid_reactions.items():
        for lipid in value['reactants']:
            acyl_chains_in_lipid = lipid.name.split(value['splitter'])[0]
            acyl_chains.add(acyl_chains_in_lipid)
    print(acyl_chains)







if __name__ == '__main__':
    model = load()
    e_lipid_reactions = {"DAG": {"splitter": "-sn-glycerol", "reactants": model.reactions.get_by_id("e_DAG__in").reactants}}
    filter_components_by_biomass(e_lipid_reactions)

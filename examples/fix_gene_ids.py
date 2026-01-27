
import sys
from cobra.io import read_sbml_model, write_sbml_model


def fix_gene_ids(model):
    """
    Fix gene ids in the model
    Parameters
    ----------
    model: MyModel
        Model to fix

    Returns
    -------
    model: MyModel
        Model with fixed gene ids
    """
    for gene in model.genes:
        if gene.id.endswith("_"):
            gene.id = gene.id[:-1]
        if gene.name.endswith("_"):
            gene.name = gene.name[:-1]
        for reaction in gene.reactions:
            reaction.gene_reaction_rule = reaction.gene_reaction_rule.replace(gene.id + "_", gene.id)
    return model


if __name__ == '__main__':
    args = sys.argv
    model = read_sbml_model(args[1])
    new_modelpath = args[1].replace(".xml", "_fixed.xml").replace(".sbml", "_fixed.sbml")
    fixed_model = fix_gene_ids(model)
    write_sbml_model(fixed_model, new_modelpath)

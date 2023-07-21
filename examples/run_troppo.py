import math
import sys

import cobra.io

sys.path.insert(0, r"C:\Users\Bisbii\PythonProjects\ExpAlgae\src")
# sys.path.insert(0, r"/home/algae/GSMMutils/src")
from GSMMutils.model.COBRAmodel import *
import seaborn as sns
from GSMMutils.omics.troppo import integration_pipeline
from cobra.flux_analysis import find_blocked_reactions
from GSMMutils.io import read_csv
from GSMMutils.omics.omics_integration import OmicsIntegration
sns.set(rc={'figure.figsize':(35, 8.27)})



def read_model(data_directory):
    model_to_load = MyModel(join(data_directory, "models/model_with_trials.xml"), "e_Biomass__cytop")
    model_to_load.add_medium(join(data_directory, "media.xlsx"), "base_medium")
    model_to_load.exchanges.EX_C00009__dra.bounds = (-0.05, 1000)
    model_to_load.exchanges.EX_C00244__dra.bounds = (-5, 1000)
    return model_to_load


def main():
    # model = read_model("../data")
    model = MyModel(join("../data", "models/model_with_trials.xml"), "e_Biomass__cytop")
    model.add_medium(join("../data", "media.xlsx"), "media_with_starch")
    print(model.optimize())
    blocked = find_blocked_reactions(model)
    model.remove_reactions(blocked)
    # model_consistent = fastcc(model)
    # cobra.io.write_sbml_model(model_consistent, r"C:\Users\Bisbii\PythonProjects\GSMMutils\data\models\model_consistent.xml")
    # model = MyModel(r"C:\Users\Bisbii\PythonProjects\ExpAlgae\data\models\model_consistent.xml", "e_Biomass__cytop")
    omics = OmicsIntegration(r"../data/omics/output.txt", samples_names={"SRR7984026Aligned.out.sam":"LL_1",
                                                                    "SRR7984027Aligned.out.sam": "LL_2",
                                                                    "SRR7984028Aligned.out.sam":"LL_3",
                                                                    "SRR7984029Aligned.out.sam":"ML_1",
                                                                    "SRR7984030Aligned.out.sam":"ML_2",
                                                                    "SRR7984031Aligned.out.sam":"ML_3",
                                                                    "SRR7984032Aligned.out.sam":"HL_1",
                                                                    "SRR7984033Aligned.out.sam":"HL_2",
                                                                    "SRR7984034Aligned.out.sam":"HL_3",
                                                                  }, model=model)

    omics.getmm = read_csv(r"../data/omics/getmm_light.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics.sum_tech_reps()
    #
    dataset = omics.counts.applymap(lambda x: math.log2(x + 1)).T
    # dataset = omics.counts.T
    ids = [ 111,  279,  282,  311,  319,  591, 1122, 1125, 1614 ,1885, 1893, 2058, 2201, 2268, 2622, 2640, 2690, 2750, 2768, 2877, 2910 ,3175, 3305]
    to_ignore = []
    print("#" * 100)
    # to_ignore = [111, 310, 1117, 1285, 1969, 2191, 2608, 2751]
    # for reac_id in to_ignore:
    #     print(reac_id)
    #     print(model.reactions[reac_id].id, model.reactions[reac_id].name)
    reaction_ids = []
    for i in ids:
        if i not in to_ignore:
            reaction = model.reactions[i]
            reaction_ids.append(reaction.id)
    # reaction_ids += ["BMGR5416__lip", "BMGR5845__lip", "BMGR5422__lip", "DGDGS_183_164__chlo"]
    reaction_ids += ["DGDGS_183_164__chlo"]
    # reaction_ids = []
    to_remove = []
    with model:
        i = 0
        for reaction_id in reaction_ids:
            if i< len(ids):
                print(ids[i], reaction_id, model.reactions.get_by_id(reaction_id).name)
            else:
                print(reaction_id, model.reactions.get_by_id(reaction_id).name)
            model.remove_reactions([reaction_id])
            to_remove.append(reaction_id)
            i+=1
    # model.remove_reactions(to_remove)
    # blocked = find_blocked_reactions(model)
    # print(blocked)
    # print(len(blocked))
    # model.remove_reactions(blocked)
    # print(model.reactions[975]) #[ 126 1867 2211]
    # print(model.reactions[3093])
    # print(model.reactions[1867])
    # print(model.reactions[3070])
    # model.write(r"C:\Users\Bisbii\PythonProjects\ExpAlgae\data\models\consistent_model.xml")
    # print(model.optimize())
    # print(dataset.quantile(0.75, axis=1))
    integration_pipeline(dataset, r"..\data\omics\light", "fastcore", 2.8, 8, model.model, output_dir="../data/omics/light")


def main_reduced_model():
    model = MyModel(join("../data", "models/model_with_trials.xml"), "e_Biomass__cytop")
    model.add_medium(join("../data", "media.xlsx"), "media_with_starch")
    consistent_model = cobra.flux_analysis.fastcc(model, flux_threshold=2)
    omics = OmicsIntegration(r"../data/omics/output.txt", samples_names={"SRR7984026Aligned.out.sam":"LL_1",
                                                                        "SRR7984027Aligned.out.sam": "LL_2",
                                                                        "SRR7984028Aligned.out.sam":"LL_3",
                                                                        "SRR7984029Aligned.out.sam":"ML_1",
                                                                        "SRR7984030Aligned.out.sam":"ML_2",
                                                                        "SRR7984031Aligned.out.sam":"ML_3",
                                                                        "SRR7984032Aligned.out.sam":"HL_1",
                                                                        "SRR7984033Aligned.out.sam":"HL_2",
                                                                        "SRR7984034Aligned.out.sam":"HL_3",
                                                                  }, model=model)

    omics.getmm = read_csv(r"../data/omics/getmm_light.tsv", index_name='GeneID', index_col=0, comment='#', sep='\t')
    omics.sum_tech_reps()
    dataset = omics.counts.applymap(lambda x: math.log2(x + 1)).T
    print(consistent_model.optimize())
    print(dataset.quantile(0.75, axis=1))
    print(len(consistent_model.reactions))
    blocked = find_blocked_reactions(consistent_model)
    print(len(blocked))
    integration_pipeline(dataset, r"..\data\omics\light", "fastcore", 2.3, 8, consistent_model, output_dir="../data/omics/light")

def get_non_constant_columns(df):
    constant_columns = []
    for col in df.columns:
        if len(df[col].unique()) != 1:
            constant_columns.append(col)
    return constant_columns
def get_different_reactions():
    model = MyModel(join("../data", "models/model_consistent.xml"), "e_Biomass__cytop")
    data = pd.read_csv(r"..\data\omics\light\results.csv", index_col=0)
    constant_columns = get_non_constant_columns(data)
    model.get_reactions_pathways_map()
    pathway_counts = {}
    for column in constant_columns:
        for pathway in model.reactions_pathways_map[column]:
            if pathway not in pathway_counts:
                pathway_counts[pathway] = 1
            else:
                pathway_counts[pathway] += 1
    print(pathway_counts)




if __name__ == '__main__':
    # model = MyModel(join("../data", "models/model.xml"), "e_Biomass__cytop")
    # model.add_medium(join("../data", "media.xlsx"), "media_with_starch")
    # model.write(r"C:\Users\Bisbii\PythonProjects\GSMMutils\data\models\model_with_media.xml")
    main()
    # get_different_reactions()
    # main_reduced_model()
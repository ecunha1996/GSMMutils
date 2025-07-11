import glob
from collections import defaultdict

import pandas as pd
from Bio import SeqIO

LOCTREE_DEEPLOC_MAP = {"cytoplasm": "Cytoplasm", "secreted": "Extracellular", "mitochondrion": "Mitochondrion", "nucleus": "Nucleus", "peroxisome": "Peroxisome",
    "chloroplast": "Mitochondrion", "endoplasmic reticulum": "Endoplasmic reticulum", "lysosome": "Lysosome", "golgi apparatus": "Golgi apparatus", "vacuole": "Lysosome/Vacuole",
                       "plastid": "Mitochondrion", "plasma membrane": "Cell membrane", "endoplasmic reticulum membrane": "Endoplasmic reticulum", "mitochondrion membrane": "Mitochondrion",
                       "peroxisome membrane": "Peroxisome", "golgi apparatus membrane": "Golgi apparatus", "lysosome membrane": "Lysosome", "vacuole membrane": "Lysosome/Vacuole",
                       "chloroplast membrane": "Mitochondrion", "nucleus membrane": "Nucleus"}



def load_deeploc(path):
    """
    Load DeepLoc annotation files from a specified path.
    Parameters
    ----------
    path: str
        Path to the DeepLoc annotation files.

    Returns
    -------
    deeploc: pd.DataFrame
        DeepLoc annotation data.
    """
    # load csv and then merge them
    deeploc = pd.concat([pd.read_csv(f, index_col=0) for f in glob.glob(path + '/*.csv')])
    deeploc = deeploc[~deeploc.index.duplicated(keep='first')]
    return deeploc

def assert_deeploc(deeploc):
    fasta_file = r"C:\Users\Bisbii\Desktop\compartments_juliana\Bjead1_1.faa"
    protein_ids = []
    handle = SeqIO.parse(fasta_file, "fasta")
    biggest_protein, biggest_protein_len = None, 0
    for record in handle:
        if record.id not in deeploc.index:
            protein_ids.append(record)
            if len(record.seq) > biggest_protein_len:
                biggest_protein = record.id
                biggest_protein_len = len(record.seq)
    print(biggest_protein, biggest_protein_len)
    print(len(protein_ids))
    protein_ids = [record for record in protein_ids if record.id != biggest_protein]
    # write in fasta
    SeqIO.write([record for record in protein_ids], r"C:\Users\Bisbii\Desktop\compartments_juliana\missing_proteins.fasta", "fasta")




def load_loctree3(path):
    """
    Load LocTree3 annotation file from a specified path.
    Parameters
    ----------
    path: str
        Path to the LocTree3 annotation file.

    Returns
    -------
    loctree3: pd.DataFrame
        LocTree3 annotation data.
    """
    return pd.read_csv(path, sep='\t').drop(['Gene Ontology Terms'], axis=1)


def get_stats_deeploc(deeploc: pd.DataFrame):
    """
    Get statistics of the DeepLoc annotation data.
    Parameters
    ----------
    deeploc

    Returns
    -------

    """
    deeploc = deeploc[~deeploc.index.duplicated(keep='first')]
    print(deeploc.shape)
    data_as_dict = deeploc.to_dict(orient='index')
    localizations_counts = {}
    for key, value in data_as_dict.items():
        if isinstance(value['Localizations'], str):
            localization = value['Localizations'].split('|')
        else:
            localization = value['Localizations']
        for loc in localization:
            if loc not in localizations_counts:
                localizations_counts[loc] = 0
            localizations_counts[loc] += 1
    print(localizations_counts)


def merge_results(deeploc, loctree):
    """
    This function will check wich proteins are predicted to be in the nucleous (bias of DeepLoc) and then will check what was the prediction made by LocTree3.
    If locTree3 predicted the protein to be in the nucleous, there will be no changes.
    If LocTree3 predicted the protein to be in another compartment, the protein will be updated in the deeploc table with a new prediction.
    The Localizations column will be updated.
    Parameters
    ----------
    deeploc
    loctree

    Returns
    -------

    """
    counter = 0
    # check if the protein is in the nucleous according to deeploc
    deeploc['Localizations'] = deeploc['Localizations'].apply(lambda x: x.split('|'))
    proteins_in_nucleous = deeploc[deeploc['Localizations'].apply(lambda x: 'Nucleus' in x)].index.tolist()
    # check if the protein is in the nucleous according to loctree
    for protein in proteins_in_nucleous:
        if protein in loctree['Protein Id'].tolist():
            loctree_prediction = loctree[loctree['Protein Id'] == protein]
            loctree_loc, score = loctree_prediction['Localization'].values[0], loctree_prediction['Score'].values[0]
            loctree_loc = LOCTREE_DEEPLOC_MAP.get(loctree_loc, loctree_loc)
            # print(f"Protein {protein} is in the nucleous according to DeepLoc, but LocTree3 predicted it to be in {loctree_loc} with a score of {score}")
            if loctree_loc != "nucleus" and loctree_loc != "Nucleus" and loctree_loc != "Extracellular":
                deeploc.at[protein, 'Localizations'].append(loctree_loc)
                deeploc.at[protein, 'Localizations'].remove('Nucleus')
                counter+=1
    # merge strings in column Localizations using | as separator
    deeploc['Localizations'] = deeploc['Localizations'].apply(lambda x: '|'.join(x))
    return deeploc


def get_gtm_from_interpro():
    """
    Get the Gene Tree Match (GTM) from InterPro.
    Returns
    -------

    """
    from collections import defaultdict
    df = pd.read_csv(r"C:\Users\Bisbii\PythonProjects\genome-analysis\results\interproscan\Stramenopiles\we_3730_nuc\we_3730_nuc.tsv", sep="\t", header=None)
    # Define column names (as identified earlier)
    df.columns = ["Gene", "Checksum", "Seq_Length", "Database", "Signature_Acc",
                  "Description", "Start", "End", "E-value", "Status", "Date",
                  "InterPro_Acc", "InterPro_Desc"]

    # Filter out rows without InterPro IDs
    df = df.dropna(subset=["InterPro_Acc"])

    # Create a dictionary where keys are InterPro IDs and values are lists of genes
    interpro_to_genes = defaultdict(set)
    for _, row in df.iterrows():
        if row["InterPro_Acc"] != "-":
            #add interpro acc, description, and gene to the dictionary
            # interpro_to_genes[row["InterPro_Acc"]].add(row["Gene"])
            interpro_to_genes[row["InterPro_Desc"]].add(row["Gene"])




    # Write the GMT file
    gmt_file = r"C:\Users\Bisbii\PythonProjects\omics-integration\results\ngaditana\PRJNA589063\deg\interproscan_annotations.gmt"
    with open(gmt_file, "w") as f:
        for interpro_id, genes in interpro_to_genes.items():
            # Convert gene set to a tab-separated string
            gene_list = "\t".join(genes)
            # Use InterPro ID as the pathway name and provide a placeholder description
            f.write(f"{interpro_id}\tInterPro Pathway\t{gene_list}\n")

    print(f"GMT file saved as {gmt_file}")

def convert_KO2GMT():
    file_path = "C:/Users/Bisbii/Downloads/user_ko(1).txt"  # Replace with your actual TSV file
    df = pd.read_csv(file_path, sep="\t", header=None, names=["Gene", "KO"], dtype=str) #.iloc[:20,:]

    # Convert to GMT format
    gmt_path = r"C:\Users\Bisbii\PythonProjects\omics-integration\results\dsalina\PRJNA495151\deg\custom_pathways_KO.gmt"
    df.fillna("", inplace=True)

    kos = df["KO"].unique()
    print(f"Found {len(kos)} unique KO terms")

    new_dict = {}

    chunk_size = 10
    for i in range(0, len(kos), chunk_size):
        chunk = kos[i:i + chunk_size]
        kegg_map = get_kegg_map(chunk)
        if kegg_map is not None:
            new_dict.update(kegg_map)

    df["KEGG Map"] = df["KO"].map(new_dict)

    # Create GMT structure
    gmt_dict = {}

    for _, row in df.iterrows():
        gene = row["Gene"]
        pathways = row["KEGG Map"]
        if not isinstance(pathways, float):
            for pathway in pathways:
                if pathway not in gmt_dict:
                    gmt_dict[pathway] = []
                gmt_dict[pathway].append(gene)

    with open(gmt_path, "w") as f:
        for pathway, genes in gmt_dict.items():
            gens_parsed = '\t'.join(genes)
            f.write(f"{pathway}\tDescription\t{gens_parsed}\n")

    print(f"GMT file saved as {gmt_path}")


global map_map_name
map_map_name = {}


def get_kegg_map(ko_accessions):
    """
    Get the KEGG map associated with a given KO accession.
    Parameters
    ----------
    ko_accession: str
        KO accession number.

    Returns
    -------
    kegg_map: str
        KEGG map associated with the KO accession.
    """
    import requests
    ko_as_str = "+".join(ko_accessions)
    url = "https://rest.kegg.jp/link/pathway/ko:" + ko_as_str
    response = requests.get(url)
    res = {}
    if response.status_code == 200:
        results = response.text.strip().split("\n")
        for map in results:
            ko, map = map.split("\t")
            if "map" in map:
                ko_parsed = ko.split(":")[1]
                map_parsed  = map.split(":")[1]
                if ko_parsed not in res:
                    res[ko_parsed] = []
                res[ko_parsed].append(map_parsed)
    all_maps = list(set([map for maps in res.values() for map in maps if map not in map_map_name.keys()]))

    chunk_size = 5
    all_maps_chunks = [all_maps[i:i + chunk_size] for i in range(0, len(all_maps), chunk_size)]
    for chunk in all_maps_chunks:
        map_url = "https://rest.kegg.jp/get/" + "+".join(chunk)
        response = requests.get(map_url)
        if response.status_code == 200:
            results = response.text.strip().split("\n")
            for line in results:
                if line.startswith("ENTRY"):
                    map_id = line.replace("ENTRY", "").strip().split()[0]
                if line.startswith("NAME"):
                    map_name= line.replace("NAME", "").strip()
                if line.startswith("CLASS"):
                    map_class = line.replace("CLASS", "").strip()
                    if "Metabolism" in map_class:
                        map_map_name[map_id] = map_name
    new_dict = {}
    for ko, maps in res.items():
        for map in maps:
            if map in map_map_name.keys():
                if ko not in new_dict.keys():
                    new_dict[ko] = []
                new_dict[ko].append(map_map_name[map])
    return new_dict



if __name__ == '__main__':
    # deeploc = load_deeploc(r"C:\Users\Bisbii\Desktop\compartments_juliana\deeploc")
    # # assert_deeploc(deeploc)
    # get_stats_deeploc(deeploc)
    # loctree = load_loctree3(r"C:\Users\Bisbii\Desktop\compartments_juliana\results.lc3")
    # merged = merge_results(deeploc, loctree)
    # get_stats_deeploc(merged)
    # merged.to_csv(r"C:\Users\Bisbii\Desktop\compartments_juliana\deeploc_merged.csv")
    # get_gtm_from_interpro()
    # convert_KO2GMT()
    from cobra.io import read_sbml_model
    # model = read_sbml_model(r"C:\Users\Bisbii\Downloads\carve_model_LB.xml")
    # print(model.exchanges.EX_o2_e.bounds)
    model_2 = read_sbml_model(r"C:\Users\Bisbii\PythonProjects\GSMMutils\data\models\iJN678.xml")
    metabolic_reactions = [reaction for reaction in model_2.reactions if reaction not in model_2.boundary and len(reaction.compartments) == 1]
    print([r.id for r in metabolic_reactions if not r.genes])
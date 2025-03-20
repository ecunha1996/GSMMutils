import json
from Bio.KEGG.REST import kegg_get
from KEGG_parser.parsers import parse_pathway
from requests import get


aliases = {"Glycolysis / Gluconeogenesis": "Glycolysis",
           "Carbon fixation in photosynthetic organisms": "Carbon fixation by Calvin cycle",
           "Carbon fixation pathways in prokaryotes": "Other carbon fixation pathways",
            "Glycine, serine and threonine metabolism": "Glycine serine and threonine metabolism",
           "Porphyrin and chlorophyll metabolism": "Porphyrin metabolism"
           }



def search_pathway_map_id(pathway_name):
    """
    Search for a pathway in KEGG and return the pathway id
    Parameters
    ----------
    pathway_name: str
        Name of the pathway to search for
    Returns
    -------
    pathway_id: str
        KEGG id of the pathway
    """
    url = f'https://rest.kegg.jp/find/pathway/{pathway_name.replace("/", "").replace(",", "")}'
    response = get(url)
    if response.text == "\n":
        pathway_name = aliases.get(pathway_name, pathway_name)
        url = f'https://rest.kegg.jp/find/pathway/{pathway_name.replace("/", "").replace(",", "")}'
        response = get(url)
    for line in response.text.splitlines():
        if len(line.split('\t')) == 1:
            continue
        pathway_id, pathway_name_from_kegg = line.split('\t')
        if pathway_name == pathway_name_from_kegg:
            return pathway_id.strip("path:")


def get_kegg_pathways(pathway_id=None):
    """
    Get the pathways from KEGG
    Parameters
    ----------
    pathway_id: str
        KEGG id of the pathway

    Returns
    -------
    result: dict
        Dictionary with the pathways
    """
    if pathway_id:
        try:
            result = kegg_get(pathway_id)
        except Exception as e:
            print(pathway_id)
            print(e)
            return None
        if result:
            result = result.read()
        result = parse_pathway(result.replace("///", ""))
        return result


def create_related_pathways_map(universal_model, folder_path):
    """
    Create a map of related pathways for each pathway in the universal model
    Parameters
    ----------
    universal_model: MyModel
        Universal model
    folder_path: str
        Path to the folder where the map will be saved

    Returns
    -------
    res: dict
        Dictionary with the related pathways for each pathway in the universal model
    """
    res = {}
    for pathway in universal_model.groups:
        related_pathways = get_related_pathways(pathway.name, res)
        res[pathway.name] = related_pathways
    json.dump(res, open(f'{folder_path}/../related_pathways_map.json', 'w'))
    return res


def get_related_pathways(pathway_name, related_pathways_map):
    """
    Get the related pathways for a pathway
    Parameters
    ----------
    pathway_name: str
        Name of the pathway
    related_pathways_map: dict
        Dictionary with the related pathways for each pathway in the universal model

    Returns
    -------
    pathways: list
        List of related pathways
    """
    if pathway_name in related_pathways_map:
        return related_pathways_map[pathway_name]
    pathways = []
    try:
        pathway_id = search_pathway_map_id(pathway_name)
        if pathway_id is not None:
            kegg_result = get_kegg_pathways(pathway_id)
            for pathway in kegg_result['REL_PATHWAY']:
                pathways.append(pathway[1])
    except Exception as e:
        print(e)
    return pathways

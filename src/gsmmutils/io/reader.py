import pandas as pd

def read_matrix(filename, **kwargs):
    matrix = pd.read_excel(filename, sheet_name=None, **kwargs)
    for value in matrix.values():
        value.index = value.index.astype(str)
    return matrix


def read_simulation():
    # TODO: Implement
    pass


def read_csv(filename, index_name=None, **kwargs):
    data = pd.read_csv(filename, **kwargs)
    data.index = data.index.astype(str)
    data.index = [index_name.split(".")[0] for index_name in data.index]
    data = data.rename_axis(index_name)
    return data


def read_excel(filename, index_name=None, **kwargs):
    data = pd.read_excel(filename, **kwargs, engine="openpyxl")
    for sheet in data:
        data[sheet].index = data[sheet].index.astype(str)
        data[sheet].index = [index_name.split(".")[0] for index_name in data[sheet].index]
        data[sheet] = data[sheet].rename_axis(index_name)
    return data

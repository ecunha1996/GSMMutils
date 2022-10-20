
import pandas as pd



def read_matrix(filename, **kwargs):
    matrix = pd.read_excel(filename, **kwargs)
    for value in matrix.values():
        value.index = value.index.astype(str)
    return matrix



def read_simulation():
    pass



def read_csv(filename, **kwargs):
    data = pd.read_csv(filename, **kwargs)
    data.index = data.index.astype(str)
    return data
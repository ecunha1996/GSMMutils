
import pandas as pd



def read_matrix(filename, **kwargs):
    matrix = pd.read_excel(filename, **kwargs)
    return matrix







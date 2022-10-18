
import pandas as pd



def write_matrix(matrix, filename, **kwargs):
    matrix.to_excel(filename, **kwargs)

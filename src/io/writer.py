
import pandas as pd



def write_matrix(matrix, filename, **kwargs):
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')
    for df_name, df in matrix.items():
        index_name = df.index.name
        df.insert(0, index_name, df.index)
        df.to_excel(writer, sheet_name=df_name, index=False, **kwargs)
        worksheet = writer.sheets[df_name]
        for idx, col in enumerate(df):
            series = df[col]
            max_len = max((
                series.astype(str).map(len).max(),
                len(str(series.name))
            )) + 3
            worksheet.set_column(idx, idx, max_len)
    writer.save()


def write_simulation(results, filename, **kwargs):
    writer = pd.ExcelWriter(f"{filename}.xlsx", engine='xlsxwriter')
    for df_name, df in results.items():
        df = pd.DataFrame(df.fluxes)
        index_name = df.index.name or "reaction_id"
        df.insert(0, index_name, df.index)
        df.to_excel(writer, sheet_name=df_name, index=False, **kwargs)
        worksheet = writer.sheets[df_name]
        for idx, col in enumerate(df):
            series = df[col]
            max_len = max((
                series.astype(str).map(len).max(),
                len(str(series.name))
            )) + 3
            worksheet.set_column(idx, idx, max_len)
    writer.save()
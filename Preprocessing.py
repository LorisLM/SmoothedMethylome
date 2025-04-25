import pandas as pd


def merge(df1, df2):
    redundant_columns = df1.columns.intersection(df2.columns)

    df2 = df2.drop(columns=redundant_columns, errors='ignore')
    df2.rename(columns={'Markers': 'Predictor'}, inplace=True)

    df_merged = pd.merge(df1, df2, on='Predictor', how='left')
    df_merged['Chromosome'] = df_merged['Chromosome'].astype(str).str.replace('.0', '', regex=False).str.strip()

    return df_merged


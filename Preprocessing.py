import pandas as pd


def merge(df1, df2):
    redundant_columns = df1.columns.intersection(df2.columns)

    df2 = df2.drop(columns=redundant_columns, errors='ignore')
    df2.rename(columns={'Markers': 'Predictor'}, inplace=True)

    df_merged = pd.merge(df1, df2, on='Predictor', how='left')

    return df_merged


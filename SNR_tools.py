import numpy as np


def signal_noise_gene(df, gene, target_method):

    df_gene = df[df['Gene Name'] == gene]
    signal = df_gene[target_method].mean()

    df_filtered = df[(df['Gene Name'] != gene) | (df['Gene Name'].isnull())]
    noise = df_filtered[target_method].mean()

    return signal/noise

def signal_noise_bootstrap_gene(df, gene, target_method, noise_percentage=0.05):

    df_gene = df[df['Gene Name'] == gene]
    signal = df_gene[target_method].mean()

    df_filtered = df[df['Gene Name'] != gene]
    noise_values = df_filtered[target_method].dropna().values

    global_noise_sample = np.random.choice(noise_values, size=int(noise_percentage * len(noise_values)), replace=True)
    noise = np.mean(global_noise_sample)

    return signal / noise


def signal_noise_pangenomic(df, window_size=3, target_method='T_log10_P_smoothed'):
    df_sorted = df.sort_values(by=['Chromosome', 'Position']).copy()

    df_sorted['Signal'] = df_sorted[target_method].rolling(window=window_size, center=True, min_periods=1).mean()

    global_noise = df_sorted[target_method].mean()

    df_sorted['Signal_Noise_Ratio'] = df_sorted['Signal'] / global_noise

    return df_sorted[['Predictor', 'Chromosome', 'Position', 'Gene Name', 'Signal_Noise_Ratio']].sort_values(by="Signal_Noise_Ratio", ascending=False)
import numpy as np


def signal_noise_for_gene(df, gene, smooth_column, bootstrap=False, noise_percentage=0.05):
    df_gene = df[df['Gene Name'] == gene]
    signal = df_gene[smooth_column].mean()

    df_noise = df[df['Gene Name'] != gene]

    if bootstrap:
        noise_values = df_noise[smooth_column].dropna().values

        sample_size = max(1, int(noise_percentage * len(noise_values)))
        noise_sample = np.random.choice(noise_values, size=sample_size, replace=True)
        noise = noise_sample.mean()
    else:
        df_noise = df_noise.copy()
        noise = df_noise[smooth_column].mean()

    return signal / noise


def signal_noise_pangenomic(df, smooth_column,window_size=3):
    df_sorted = df.sort_values(by=['Chromosome', 'Position']).copy()

    df_sorted['Signal'] = df_sorted[smooth_column].rolling(window=window_size, center=True, min_periods=1).mean()

    global_noise = df_sorted[smooth_column].mean()

    df_sorted['Signal_Noise_Ratio'] = df_sorted['Signal'] / global_noise

    return df_sorted
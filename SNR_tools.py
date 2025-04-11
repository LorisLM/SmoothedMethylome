import numpy as np


def signal_noise_gene(df, gene, target_method='T_log10_P_smoothed'):
    """
      Calculates a simple signal-to-noise ratio for a specific gene.

      The signal is defined as the mean value of the 'target_method' for the given 'gene'.
      The noise is defined as the mean value of the 'target_method' for all other genes.

      Args:
          df (pd.DataFrame): DataFrame with 'Gene Name' and the target method column.
          gene (str): The name of the gene to consider as signal.
          target_method (str, optional): The column name containing the values to analyze.
                                         Defaults to 'T_log10_P_smoothed'.

      Returns:
          float: The signal-to-noise ratio for the specified gene.
      """

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
    """
    This function computes a genome-wide signal-to-noise ratio.
    It smooths a chosen metric across the genome using a rolling mean and then compares this smoothed signal to the overall average of that metric.

    Args:
        df (pd.DataFrame): DataFrame containing genomic data with 'Chromosome', 'Position',
                           'Gene Name', 'Predictor', and the 'target_method' column.
        window_size (int, optional): The size of the rolling window for smoothing.
                                     Must be a positive integer. Defaults to 3.
        target_method (str, optional): The name of the column to use for signal calculation.
                                       Defaults to 'T_log10_P_smoothed'.

    Returns:
        pd.DataFrame: DataFrame containing 'Predictor', 'Chromosome', 'Position',
                      'Gene Name', and the calculated 'Signal_Noise_Ratio', sorted
                      in descending order of 'Signal_Noise_Ratio'.
    """
    df_sorted = df.sort_values(by=['Chromosome', 'Position']).copy()

    df_sorted['Signal'] = df_sorted[target_method].rolling(window=window_size, center=True, min_periods=1).mean()

    global_noise = df_sorted[target_method].mean()

    df_sorted['Signal_Noise_Ratio'] = df_sorted['Signal'] / global_noise

    return df_sorted[['Predictor', 'Chromosome', 'Position', 'Gene Name', 'Signal_Noise_Ratio']].sort_values(by="Signal_Noise_Ratio", ascending=False)
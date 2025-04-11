import numpy as np
import pandas as pd

from SNR_tools import signal_noise_bootstrap_gene
from Smoothing import mean_smoothing, savgol_smoothing


def bootstrap_distribution(df):
    """
    Generates a DataFrame containing bootstrapped signal-to-noise ratio distributions
    for the gene 'MMACHC' across different mean smoothing window sizes.

    This function iterates through odd window sizes from 1 to 19. For each window size,
    it performs 1000 bootstrap resamplings of the input DataFrame. In each bootstrap
    iteration, it calculates a signal-to-noise ratio for the gene 'MMACHC' using the
    'signal_noise_bootstrap_gene' function, with 'T_log10_P_smoothed' as the target method.

    Args:
        df (pd.DataFrame): DataFrame with 'Gene Name' and 'T_log10_P' columns.

    Returns:
        pd.DataFrame: DataFrame with columns 'n' and 'signal_noise_bootstraps'.
                     - 'n' (int): The mean smoothing window size (1, 3, 5, ..., 19).
                     - 'signal_noise_bootstraps' (list of float): A list of 1000
                       bootstrapped signal-to-noise ratios for 'MMACHC' for the
                       corresponding window size.
    """
    results = []
    smoothed_versions = {n: mean_smoothing(df[['Gene Name', 'T_log10_P']].copy(), n) for n in range(1, 20, 2)}

    for n, smoothed_df in smoothed_versions.items():
        signal_noise_bootstraps = []

        for _ in range(1000):
            ratio = signal_noise_bootstrap_gene(smoothed_df, gene='MMACHC', target_method='T_log10_P_smoothed')
            signal_noise_bootstraps.append(ratio)

        results.append({'n': n,'signal_noise_bootstraps': signal_noise_bootstraps})

    return pd.DataFrame(results)

def bootstrap_mean_with_CI(df, gene):
    """
    Performs signal-to-noise ratio estimation via bootstrapping for different smoothing methods
    (mean and Savitzky-Golay) and parameters, then returns the results in a DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing gene information.
                           It should have at least a column named 'Gene Name'
                           and the 'T_log10_P'.
        gene (str): The name of the gene for which to calculate the bootstrapped
                    mean and confidence interval.


    Returns:
        pd.DataFrame: A dataframe containing the bootstrapped mean ('mean'), the lower
                      bound of the confidence interval ('lower_ci'), and the upper bound
                      of the confidence interval ('upper_ci') for the specified gene.
    """
    results = []

    # Mean smoothing
    for n in range(1, 20, 2):
        smoothed_df = mean_smoothing(df[['Gene Name', 'T_log10_P']].copy(), n)
        mean_bootstrap, conf_interval = bootstrap_with_confidence(smoothed_df, gene, 'T_log10_P_smoothed')

        results.append({
            'n': n,
            'order': 0,
            'method': 'mean',
            'signal_noise_mean': mean_bootstrap,
            'conf_lower': conf_interval[0],
            'conf_upper': conf_interval[1]
        })

    # Savitzky-Golay smoothing (orders 1 and 2)
    for order in [1, 2, 4]:
        for n in range(3, 20, 2):
            if n > order:
                df_copy = df.copy()
                savgol_smoothing(df_copy, n, order)
                mean_bootstrap, conf_interval = bootstrap_with_confidence(df_copy, gene, 'T_log10_P_smoothed')

                results.append({
                    'n': n,
                    'order': order,
                    'method': 'savgol',
                    'signal_noise_mean': mean_bootstrap,
                    'conf_lower': conf_interval[0],
                    'conf_upper': conf_interval[1]
                })
    pd.DataFrame(results).to_excel("smooth_CI.xlsx")
    return pd.DataFrame(results)

def bootstrap_with_confidence(df, gene, target_method, noise_percentage=0.10, n_bootstrap=750):

    bootstrap_samples = [
        signal_noise_bootstrap_gene(df, gene, target_method, noise_percentage) for _ in range(n_bootstrap)
    ]

    mean_value = np.mean(bootstrap_samples)
    conf_interval = np.percentile(bootstrap_samples, [2.5, 97.5])

    return mean_value, conf_interval
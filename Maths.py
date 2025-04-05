import numpy as np
import pandas as pd

from SNR_tools import signal_noise_bootstrap_gene
from Smoothing import mean_smoothing


def bootstrap_distribution(df):
    results = []
    smoothed_versions = {n: mean_smoothing(df[['Gene Name', 'T_log10_P']].copy(), n) for n in range(1, 20, 2)}

    for n, smoothed_df in smoothed_versions.items():
        signal_noise_bootstraps = []

        for _ in range(1000):
            ratio = signal_noise_bootstrap_gene(smoothed_df, gene='MMACHC', target_method='T_log10_P_smoothed')
            signal_noise_bootstraps.append(ratio)

        results.append({'n': n,'signal_noise_bootstraps': signal_noise_bootstraps})

    return pd.DataFrame(results)

def bootstrap_mean_with_CI(df, gene, target_method):

    results = []
    smoothed_versions = {
        n: mean_smoothing(df[['Gene Name', 'T_log10_P']].copy(), n)
        for n in range(1, 20, 2)
    }

    for n, smoothed_df in smoothed_versions.items():
        mean_bootstrap, conf_interval = bootstrap_with_confidence(smoothed_df, gene, 'T_log10_P_smoothed')

        results.append({
            'n': n,
            'signal_noise_mean': mean_bootstrap,
            'conf_lower': conf_interval[0],
            'conf_upper': conf_interval[1]
        })

    return pd.DataFrame(results)

def bootstrap_with_confidence(df, gene, target_method, noise_percentage=0.10, n_bootstrap=1000):

    bootstrap_samples = [
        signal_noise_bootstrap_gene(df, gene, target_method, noise_percentage) for _ in range(n_bootstrap)
    ]

    mean_value = np.mean(bootstrap_samples)
    conf_interval = np.percentile(bootstrap_samples, [2.5, 97.5])

    return mean_value, conf_interval
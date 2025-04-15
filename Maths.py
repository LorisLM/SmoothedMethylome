import numpy as np
import pandas as pd

from SNR_tools import signal_noise_for_gene
from Smoothing import mean_smoothing, savgol_smoothing


def bootstrap_distribution(df, gene, smooth_column, windows_size=20, n_iterations=1000):
    results = []

    for n in range(1, windows_size, 2):
        smoothed_df = mean_smoothing(df.copy(), n)

        signal_noise_bootstraps = []
        for _ in range(n_iterations):
            ratio = signal_noise_for_gene(smoothed_df, gene=gene, bootstrap=True)
            signal_noise_bootstraps.append(ratio)

        results.append({
            'n': n,
            'signal_noise_bootstraps': signal_noise_bootstraps
        })

    return pd.DataFrame(results)

def bootstrap_mean_with_CI(df, gene, target_column, windows_size=20):
    results = []

    smoothing_methods = [
        {
            'name': 'mean',
            'orders': [0],
            'window_range': range(1, windows_size, 2),
            'apply': lambda df_copy,target_column, n, order: mean_smoothing(df_copy.copy(),target_column, n)
        },
        {
            'name': 'savgol',
            'orders': [1, 2, 4],
            'window_range': range(3, windows_size, 2),
            'apply': lambda df_copy,target_column, n, order: savgol_smoothing(df_copy.copy(), target_column, n, order)
        }
    ]

    for method in smoothing_methods:
        for order in method['orders']:
            for n in method['window_range']:
                if method['name'] == 'savgol' and n <= order:
                    continue
                smoothed_df = method['apply'](df, target_column, n, order)

                mean_bootstrap, conf_interval = bootstrap_with_confidence(
                    smoothed_df, gene, 'smooth_result'
                )

                results.append({
                    'n': n,
                    'order': order,
                    'method': method['name'],
                    'signal_noise_mean': mean_bootstrap,
                    'conf_lower': conf_interval[0],
                    'conf_upper': conf_interval[1]
                })

    result_df = pd.DataFrame(results)


    return result_df

def bootstrap_with_confidence(df, gene, smooth_column, noise_percentage=0.10, n_bootstrap=20):

    bootstrap_samples = [
        signal_noise_for_gene(df, gene, smooth_column, True, noise_percentage) for _ in range(n_bootstrap)
    ]

    mean_value = np.mean(bootstrap_samples)
    conf_interval = np.percentile(bootstrap_samples, [2.5, 97.5])

    return mean_value, conf_interval
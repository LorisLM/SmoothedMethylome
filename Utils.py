import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from SNR_tools import signal_noise_gene
from Smoothing import savgol_smoothing


def savgol_looker(df_merged):
    results = []

    for n in range(3, 101, 2):
        for order in range(1, min(n, 10)):
            savgol_smoothing(df_merged, n, order)
            signal_noise_ratio = signal_noise_gene(df_merged, 'MMACHC', 'T_log10_P_smoothed')
            results.append({'n': n, 'order': order, 'signal_noise': signal_noise_ratio})

    df_results = pd.DataFrame(results)

    plt.figure(figsize=(12, 8))

    for order_val in df_results['order'].unique():
        subset = df_results[df_results['order'] == order_val]
        sns.lineplot(x='n', y='signal_noise', data=subset, label=f'order={order_val}', marker='o')

    plt.title("Signal-to-Noise Ratio vs. Window Size (n) for Different Orders")
    plt.xlabel("Window Size (n)")
    plt.ylabel("Signal-to-Noise Ratio")
    plt.legend(title="Order")

    odd_ticks = np.arange(min(df_results['n']), max(df_results['n']) + 1, 2)  # Nombres impairs
    plt.xticks(odd_ticks)
    plt.grid(axis='x', which='major')
    plt.grid(axis='y')

    plt.show()

def chromosome_sort(chrom):
    if chrom.isdigit():
        return int(chrom)
    elif chrom == "X":
        return 23
    elif chrom == "Y":
        return 24
    else:
        return 25


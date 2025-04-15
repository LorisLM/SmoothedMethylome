import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from SNR_tools import signal_noise_for_gene
from Smoothing import savgol_smoothing


def savgol_looker(df_merged):
    plt.rcParams['font.family'] = 'Arial'

    results = []

    for n in range(3, 101, 2):
        for order in range(1, min(n, 10)):
            savgol_smoothing(df_merged, n, order)
            signal_noise_ratio = signal_noise_for_gene(df_merged, 'MMACHC', 'smooth_result')
            results.append({'n': n, 'order': order, 'signal_noise': signal_noise_ratio})

    df_results = pd.DataFrame(results)

    plt.figure(figsize=(12, 6))
    for order_val in df_results['order'].unique():
        subset = df_results[df_results['order'] == order_val]
        sns.lineplot(x='n', y='signal_noise', data=subset, label=f'order={order_val}', marker='o')

    plt.xlabel("Windows size (width*2 + 1)", fontsize=14, labelpad=17)
    plt.ylabel("Signal-to-Noise Ratio", fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=14)

    plt.legend(
        title="Order",
        bbox_to_anchor=(1, 0.90),
        loc='upper left',
        borderaxespad=0.
    )

    sns.despine(top=True, right=True)

    odd_ticks = np.arange(min(df_results['n']), max(df_results['n']) + 1, 2)
    plt.xticks(odd_ticks)

    plt.tight_layout()

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


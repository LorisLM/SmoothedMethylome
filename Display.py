import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm

from Utils import chromosome_sort
from sklearn.cluster import DBSCAN


def plot_manhattan(df, smooth_column):
    """Generates a Manhattan plot from a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame containing the genomic data.

    Returns:
        None: Displays the Manhattan plot using matplotlib.pyplot.
    """

    df = df.sort_values(["Chromosome", "Position"])
    df["Chromosome"] = df["Chromosome"].astype(str)
    df["Position"] = df["Position"].astype(int)

    chrom_offsets = {}
    current_offset = 0
    sorted_chromosomes = sorted(df["Chromosome"].unique(), key=chromosome_sort)

    for chrom in sorted_chromosomes:
        chrom_offsets[chrom] = current_offset
        max_pos = df[df["Chromosome"] == chrom]["Position"].max()
        current_offset += max_pos

    df["Pos_cum"] = df.apply(lambda row: row["Position"] + chrom_offsets[row["Chromosome"]], axis=1)

    num_colors = len(df["Chromosome"].unique())
    colors = plt.cm.plasma(np.linspace(0, 1, num_colors))

    chrom_color_map = {chrom: colors[i] for i, chrom in enumerate(sorted_chromosomes)}
    df["Color"] = df["Chromosome"].map(chrom_color_map)

    plt.figure(figsize=(14, 6))
    plt.scatter(df["Pos_cum"], df[smooth_column], c=df["Color"], edgecolor=df["Color"], alpha=0.5)

    for chrom in chrom_offsets.values():
        plt.axvline(x=chrom, color='grey', linestyle='--', alpha=0.5)

    chrom_ticks = {chrom_offsets[chrom] + (df[df["Chromosome"] == chrom]["Pos_cum"].max() - chrom_offsets[chrom]) / 2: chrom for chrom in chrom_offsets.keys()}
    plt.xticks(list(chrom_ticks.keys()), list(chrom_ticks.values()))

    plt.xlabel("Position on genome")
    plt.ylabel("-log10(p-value)")

    plt.show()

def plot_bootstrap_mean_with_CI(df):
    """
    Plots the bootstrapped mean and confidence interval(s) generated by the
    bootstrap_mean_with_CI function.

    Args:
        bootstrap_results (pd.DataFrame): The output from the bootstrap_mean_with_CI function.

    Returns:
        None: Displays the plot using matplotlib.pyplot and seaborn.

    Example:
        >>> df = pd.DataFrame(data)
        >>> results = bootstrap_mean_with_CI(df)
        >>> plot_bootstrap_mean_with_CI(results)
    """
    plt.rcParams['font.family'] = 'Arial'

    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(10, 8))

    for (method, order), subset in df.groupby(['method', 'order']):
        label = f"{method} (order={int(order)})" if pd.notna(order) else method
        sns.lineplot(
            x="n", y="signal_noise_mean", data=subset,
            marker="o", label=label
        )
        ax.fill_between(
            subset['n'], subset['conf_lower'], subset['conf_upper'],
            alpha=0.2
        )

    ax.set_xticks(sorted(df['n'].unique()))
    sns.despine(top=True, right=True, left=False, bottom=False)
    plt.legend(
        title="Smoothing method",
        bbox_to_anchor=(1, 0.90),
        loc='upper left',
        borderaxespad=0.,
        fontsize=17
    )
    plt.xticks(fontsize=17)
    plt.yticks(fontsize=17)

    plt.xlabel("Window size (width*2 + 1)", fontsize=17, labelpad=17)
    plt.ylabel("Signal-to-Noise Ratio", fontsize=17)

    plt.show()

def plot_bootstrap_violon(df):
    """
    Displays a violin plot of signal-to-noise ratio distributions obtained from bootstrapping,
    across different smoothing window sizes.

    Args:
        df (pd.DataFrame): The output from the bootstrap_distribution function.

    Returns:
        None: The function directly displays a matplotlib figure.
    """

    plt.rcParams['font.family'] = 'Arial'

    violin_df = df.explode('signal_noise_bootstraps').rename(columns={'signal_noise_bootstraps': 'signal_noise'})

    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(10, 8))

    sns.violinplot(x="n", y="signal_noise", data=violin_df, width=0.8, palette="viridis", inner="quartile")

    sns.despine(top=True, right=True, left=False, bottom=False)

    plt.xlabel("Windows size (width*2 + 1)", fontsize=14, labelpad=17)
    plt.ylabel("Signal-to-Noise Ratio", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)


    plt.show()

def plot_manhattan2(df, snr_filtered):
    """
    Generates a modified Manhattan plot highlighting points based on SNR.

    The 'snr_filtered' DataFrame should contain a 'Predictor' column with the identifiers
    of the points that should be highlighted (e.g., those with high Signal-to-Noise Ratio).

    Args:
        df (pd.DataFrame): DataFrame containing the genomic data.
        snr_filtered (pd.DataFrame): DataFrame containing 'Predictor' values to highlight.
        n (int): The total number of data points (used in the plot title).

    Returns:
        None: Displays the modified Manhattan plot using matplotlib.pyplot.

    Example:
        >>> s = savgol_smoothing(df_merged, 17,2)
        >>> snr_df = signal_noise_pangenomic(s)
        >>> snr_filtered = snr_df[snr_df["Signal_Noise_Ratio"] > 5.9]
        >>> plot_manhattan2(s,snr_filtered)
    """

    plt.rcParams['font.family'] = 'Arial'

    df = df.sort_values(["Chromosome", "Position"])
    df["Chromosome"] = df["Chromosome"].astype(str)
    df["Position"] = df["Position"].astype(int)

    # Position cumulative
    chrom_offsets = {}
    current_offset = 0
    sorted_chromosomes = sorted(df["Chromosome"].unique(), key=chromosome_sort)

    for chrom in sorted_chromosomes:
        chrom_offsets[chrom] = current_offset
        max_pos = df[df["Chromosome"] == chrom]["Position"].max()
        current_offset += max_pos

    df["Pos_cum"] = df.apply(lambda row: row["Position"] + chrom_offsets[row["Chromosome"]], axis=1)

    df["is_high_snr"] = df["Predictor"].isin(snr_filtered["Predictor"])

    plt.figure(figsize=(14, 6))

    # Points normaux (gris)
    plt.scatter(
        df.loc[~df["is_high_snr"], "Pos_cum"],
        df.loc[~df["is_high_snr"], "smooth_result"],
        color="lightgray",
        alpha=0.5,
        label="Noise"
    )

    high_snr_df = df.loc[df["is_high_snr"]].copy()
    gene_series = high_snr_df["Gene Name"].astype(str)
    unique_genes = sorted(gene_series[gene_series != 'nan'].unique().tolist())
    all_genes_label = ", ".join(unique_genes)

    if not high_snr_df.empty:
        first_pos = high_snr_df["Position"].min()
        last_pos = high_snr_df["Position"].max()
        position_label = f"\nPosition: {first_pos}-{last_pos}"
    else:
        position_label = ""

    plt.scatter(
        high_snr_df["Pos_cum"],
        high_snr_df["smooth_result"],
        color="crimson",
        alpha=0.9,
        label=f"{all_genes_label}{position_label}"
    )

    for chrom in chrom_offsets.values():
        plt.axvline(x=chrom, color='grey', linestyle='--', alpha=0.3)

    # Ticks pour les chromosomes
    chrom_ticks = {
        chrom_offsets[chrom] + (df[df["Chromosome"] == chrom]["Pos_cum"].max() - chrom_offsets[chrom]) / 2: chrom
        for chrom in chrom_offsets
    }
    plt.yticks(fontsize=17)
    plt.xticks(list(chrom_ticks.keys()), list(chrom_ticks.values()), fontsize=17)



    plt.xlabel("Chromosomes", fontsize=19, labelpad=17)
    plt.ylabel("-log10(p-value)", fontsize=19)
    sns.despine(top=True, right=True)

    plt.legend(
        bbox_to_anchor=(0, 1.1),
        loc='lower left',
        borderaxespad=0.,
        fontsize = 17
    )

    plt.tight_layout()
    plt.show()


def plot_volcano(
    stats: pd.DataFrame,
    *,
    effect_col: str = "delta",
    p_col: str = "pvalue",
    effect_thresh: float | None = None,
    p_thresh: float | None = None,
    title: str | None = None,
    ax: plt.Axes | None = None,
    out_file: str | None = None,
    show: bool = True,
) -> plt.Axes:
    """Create a volcano plot (|effect| vs –log10 p).

    Parameters
    ----------
    stats
        DataFrame returned by `compute_group_stats`.
    effect_col, p_col
        Column names for the effect size (delta) and p‑value.
    effect_thresh, p_thresh
        Optional thresholds drawn as vertical / horizontal lines.
    title
        Title of the plot.
    ax
        Matplotlib axes; created if *None*.
    out_file
        If given, saves the figure (PNG, PDF…).
    show
        Call ``plt.show()`` (ignored if running in a non‑interactive backend).
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.figure

    x = stats[effect_col].astype(float)
    p = stats[p_col].astype(float)
    y = -np.log10(p)

    ax.scatter(x, y, s=8, alpha=0.6)

    if effect_thresh is not None:
        ax.axvline(effect_thresh, linestyle="--", linewidth=1)
        ax.axvline(-effect_thresh, linestyle="--", linewidth=1)
    if p_thresh is not None:
        ax.axhline(-np.log10(p_thresh), linestyle="--", linewidth=1)

    ax.set_xlabel(effect_col)
    ax.set_ylabel("-log10 p-value")
    if title:
        ax.set_title(title)
    ax.grid(True, linestyle=":", linewidth=0.5)

    if out_file:
        fig.savefig(out_file, bbox_inches="tight")

    if show:
        plt.show()

    return ax


def plot_anova_volcano(stats: pd.DataFrame):
    plt.scatter(stats["delta"], np.log10(stats["bf10"]), alpha=0.5)
    plt.axhline(np.log10(3), color='orange', linestyle='--', label='BF10 = 3')
    plt.axhline(np.log10(10), color='red', linestyle='--', label='BF10 = 10')
    plt.xlabel("Delta")
    plt.ylabel("log10(BF10)")
    plt.legend()
    plt.show()

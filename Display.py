import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from Utils import chromosome_sort


def plot_manhattan(df, n):
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
    plt.scatter(df["Pos_cum"], df["T_log10_P_smoothed"], c=df["Color"], edgecolor=df["Color"], alpha=0.5)

    for chrom in chrom_offsets.values():
        plt.axvline(x=chrom, color='grey', linestyle='--', alpha=0.5)

    chrom_ticks = {chrom_offsets[chrom] + (df[df["Chromosome"] == chrom]["Pos_cum"].max() - chrom_offsets[chrom]) / 2: chrom for chrom in chrom_offsets.keys()}
    plt.xticks(list(chrom_ticks.keys()), list(chrom_ticks.values()))

    plt.xlabel("Position on genome")
    plt.ylabel("-log10(p-value)")
    plt.title("Manhattan Plot, Moving Average Smoothing (n=3)")
    #    plt.title("Manhattan Plot, Savitzky-Golay Smoothing (order=3)")

    plt.show()

def plot_bootstrap_mean_with_CI(df_results):
    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(10, 8))

    sns.lineplot(
        x="n", y="signal_noise_mean", data=df_results,
        marker="o", label="Mean Bootstrap", color="b"
    )
    ax.fill_between(
        df_results['n'], df_results['conf_lower'], df_results['conf_upper'],
        color='b', alpha=0.2, label="95% CI"
    )

    ax.set_xticks(df_results['n'])
    ax.xaxis.grid(True)
    ax.set(ylabel="Signal-to-Noise Ratio")
    sns.despine(trim=True, left=True)
    plt.title("Bootstrap Mean & Confidence Interval")
    plt.legend()
    plt.show()

def plot_bootstrap_violon(df):
    violin_df = df.explode('signal_noise_bootstraps').rename(columns={'signal_noise_bootstraps': 'signal_noise'})

    sns.set_theme(style="ticks")
    f, ax = plt.subplots(figsize=(10, 8))

    sns.violinplot(x="n", y="signal_noise", data=violin_df, width=0.8, palette="viridis", inner="quartile")

    ax.xaxis.grid(True)
    ax.set(ylabel="Signal-to-Noise Ratio")
    sns.despine(trim=True, left=True)
    plt.title("Distribution of Signal-to-Noise Ratio (Bootstrapped)")

    plt.show()

def plot_manhattan2(df, snr_filtered, n):
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

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.scatter(
        df.loc[~df["is_high_snr"], "Pos_cum"],
        df.loc[~df["is_high_snr"], "T_log10_P_smoothed"],
        color="lightgray",
        alpha=0.5,
        label="Others"
    )

    high_snr_df = df.loc[df["is_high_snr"]].copy()
    scatter = ax.scatter(
        high_snr_df["Pos_cum"],
        high_snr_df["T_log10_P_smoothed"],
        color="crimson",
        alpha=0.9,
        label="Point of interest"
    )

    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind):
        pos = scatter.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        row_data = high_snr_df.iloc[ind["ind"][0]]
        text = f"Chromosome: {row_data['Chromosome']}\nPosition: {row_data['Position']}\nGene Name: {row_data.get('Gene Name', 'N/A')}"
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(0.4)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = scatter.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

    for chrom in chrom_offsets.values():
        ax.axvline(x=chrom, color='grey', linestyle='--', alpha=0.3)

    chrom_ticks = {
        chrom_offsets[chrom] + (df[df["Chromosome"] == chrom]["Pos_cum"].max() - chrom_offsets[chrom]) / 2: chrom
        for chrom in chrom_offsets
    }
    ax.set_xticks(list(chrom_ticks.keys()))
    ax.set_xticklabels(list(chrom_ticks.values()))

    ax.set_xlabel("Position on genome")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Manhattan Plot, Mean smoothing (n=5)")
    ax.legend()
    fig.tight_layout()
    plt.show()

import clustering
from Display import *
from Maths import *
from Microarray_handler import *
from Preprocessing import merge
from SNR_tools import signal_noise_pangenomic
from Utils import *


"""
df3 = pd.read_csv("ressources/EWAS_regression_aging_450K.csv", sep=";")
df3 = merge(df3, df2)
df3 = df3.dropna(subset=['Position', 'Chromosome'])

genes = snr_filtered["Gene Name"].dropna().unique()

genes_sorted = sorted(genes)

plot_manhattan(df3, "smooth_result")
"""


df = read_sdrf("ressources/E-GEOD-30870.sdrf.txt")

newborns_df, nona_df = build_methylome_dataframes(
        "ressources/E-GEOD-30870.sdrf.txt",
        "ages",
    )

stats = run_bayesian_methylation_test(
    newborns_df, nona_df,
    name1="Newborns",
    name2="Nonagenarians", progress=True)

stats.to_csv("ressources/Anova_Age.csv")

#plot_volcano(stats, title="EWAS Newborns vs Nonagenarians",p_thresh=0.05, effect_thresh=0.10)
plot_anova_volcano(stats)




# MMACHC 2*var
"""
df1 = pd.read_csv("ressources/MMACHC1.csv", sep=";", dtype={"Chromosome": str})
df2 = pd.read_csv("ressources/Marker1.csv", sep=";")

df_merged = merge(df1, df2)
mean_smoothing(df_merged,"T_log10_P",7)
print(df_merged["smooth_result"].var())
snr_filtered = df_merged[df_merged["smooth_result"] > 2*df_merged["smooth_result"].var()]


plot_manhattan2(df_merged,snr_filtered)

"""

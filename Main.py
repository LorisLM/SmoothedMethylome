from Display import plot_manhattan, plot_manhattan2
from Preprocessing import merge
from SNR_tools import signal_noise_pangenomic
from Smoothing import mean_smoothing
from Utils import *

df = pd.read_csv("MMACHC1.csv", sep=";", dtype={"Chromosome": str})
df2 = pd.read_csv("Marker1.csv", sep=";")

df_merged = merge(df, df2)
s = mean_smoothing(df_merged, 5)

snr_df = signal_noise_pangenomic(s)
snr_filtered = snr_df[snr_df["Signal_Noise_Ratio"] > 6]



plot_manhattan2(s,snr_filtered,3)

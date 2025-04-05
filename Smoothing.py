from scipy.signal import savgol_filter


def mean_smoothing(df, n):
    df['T_log10_P_smoothed'] = df['T_log10_P'].rolling(window=n, center=True).mean().fillna(df['T_log10_P'])
    return df

def savgol_smoothing(df, n, order):
    df['T_log10_P_smoothed'] = savgol_filter(df['T_log10_P'], n, order)
    return df
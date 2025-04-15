from scipy.signal import savgol_filter


def mean_smoothing(df, target_column, n):
    df['smooth_result'] = df[target_column].rolling(window=n, center=True).mean().fillna(df[target_column])
    return df

def savgol_smoothing(df, target_column, n, order):
    df['smooth_result'] = savgol_filter(df[target_column], n, order)
    return df
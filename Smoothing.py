from scipy.signal import savgol_filter


def mean_smoothing(df, n):
    """
    Applies mean smoothing to the 'T_log10_P' column of a DataFrame.

    Args:
        df (pd.DataFrame): DataFrame with 'Gene Name' and 'T_log10_P' columns.
        n (int): The window size for mean smoothing.

    Returns:
        pd.DataFrame: DataFrame with an added 'T_log10_P_smoothed' column.
    """
    df['T_log10_P_smoothed'] = df['T_log10_P'].rolling(window=n, center=True).mean().fillna(df['T_log10_P'])
    return df

def savgol_smoothing(df, n, order):
    """
    Applies Savitzky-Golay smoothing to the 'T_log10_P'.

    The Savitzky-Golay filter smooths data by fitting successive sub-sets of adjacent
    data points with a low-degree polynomial by the method of linear least squares.

    Args:
        df (pd.DataFrame): DataFrame with 'Gene Name' and 'T_log10_P' columns.
        n (int): The length of the filter window (i.e., the number of coefficients).
                 Must be a positive odd number.
        order (int): The order of the polynomial used to fit the samples.
                   Must be less than 'n'.

    Returns:
        pd.DataFrame: DataFrame with an added 'T_log10_P_smoothed' column
                      containing the smoothed values.
    """
    df['T_log10_P_smoothed'] = savgol_filter(df['T_log10_P'], n, order)
    return df
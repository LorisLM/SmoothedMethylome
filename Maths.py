from typing import Iterable

import numpy as np
import pandas as pd

from SNR_tools import signal_noise_for_gene
from Smoothing import mean_smoothing, savgol_smoothing
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

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

def run_methylation_ttest(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    *,
    name1: str = "group1",
    name2: str = "group2",
    min_n: int = 2,
    progress: bool = True,
) -> pd.DataFrame:
    common = df1.index.intersection(df2.index)
    g1 = df1.loc[common]
    g2 = df2.loc[common]

    mean1 = g1.mean(axis=1, skipna=True)
    mean2 = g2.mean(axis=1, skipna=True)
    delta = mean2 - mean1

    iterator = common
    bar_close = lambda: None  # default no‑op
    try:
        from tqdm import tqdm  # type: ignore
        iterator = tqdm(common, desc="run_methylation_ttest", unit="CpG")
        bar_close = iterator.close
    except ModuleNotFoundError:
        print("[run_methylation_ttest] tqdm not installed; fallback to simple progress.")
        progress = "simple"
        step = max(len(common) // 10, 1)

    pvals: list[float] = []
    for i, idx in enumerate(iterator):
        v1 = g1.loc[idx].dropna().values
        v2 = g2.loc[idx].dropna().values
        if len(v1) < min_n or len(v2) < min_n:
            pvals.append(np.nan)
        else:
            _, p = ttest_ind(v1, v2, equal_var=False, nan_policy="omit")
            pvals.append(p)

        if progress == "simple" and (i + 1) % step == 0:
            print(f"  processed {i + 1}/{len(common)} CpGs…", end="\r")

    bar_close()

    if progress == "simple":
        print()  # newline after last carriage return

    out = pd.DataFrame({
        f"{name1}_mean": mean1,
        f"{name2}_mean": mean2,
        "delta": delta,
        "pvalue": pvals,
    })

    p_array = np.asarray(pvals, dtype=float)
    nan_mask = np.isnan(p_array)
    p_filled = np.where(nan_mask, 1.0, p_array)
    qvals = multipletests(p_filled, method="fdr_bh")[1]
    qvals[nan_mask] = np.nan
    out["qvalue"] = qvals

    return out

def _bayes_ttest(group1: np.ndarray, group2: np.ndarray) -> float:
    try:
        import pingouin as pg  # local import to keep dependency optional
    except ModuleNotFoundError as err:
        raise ImportError(
            "run_bayesian_methylation_test requires the 'pingouin' package.\n"
            "Install it with:  pip install pingouin"
        ) from err

    res = pg.ttest(group1, group2, correction=False)
    return float(res.loc["T-test", "BF10"])


def run_bayesian_methylation_test(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    *,
    name1: str = "Grp1",
    name2: str = "Grp2",
    progress: bool | str = False,
) -> pd.DataFrame:
    """Compute per‑CpG Bayes‑factor between two groups of β‑values.

    Parameters
    ----------
    df1, df2 : DataFrame
        CpG × samples matrices (index = CpG IDs). Must share the same dtype
        and be on the same scale (β‑values).
    name1, name2 : str
        Labels used to name the output columns.
    progress : bool or str
        If truthy, wrap the loop with `tqdm` progress bar. If a string, that
        string is used as the *desc* argument for tqdm.

    Returns
    -------
    DataFrame
        Columns: ``<name1>_mean, <name2>_mean, delta, bf10`` with CpG IDs as
        index.
    """
    if df1.index.nlevels != 1 or df2.index.nlevels != 1:
        raise ValueError("DataFrames must have CpG IDs as a flat index.")

    common = df1.index.intersection(df2.index)
    if common.empty:
        raise ValueError("df1 and df2 share no CpG IDs.")

    g1 = df1.loc[common].to_numpy(dtype=float)
    g2 = df2.loc[common].to_numpy(dtype=float)

    means1 = np.nanmean(g1, axis=1)
    means2 = np.nanmean(g2, axis=1)
    delta  = means2 - means1

    iterator: Iterable[int]
    if progress:
        try:
            from tqdm import tqdm

            desc = progress if isinstance(progress, str) else "Compute BF10"
            iterator = tqdm(range(len(common)), desc=desc)
        except ModuleNotFoundError:
            print("[run_bayesian_methylation_test] tqdm not installed – running without progress bar")
            iterator = range(len(common))
    else:
        iterator = range(len(common))

    bf10 = np.empty(len(common), dtype=float)
    for i in iterator:
        bf10[i] = _bayes_ttest(g1[i, :], g2[i, :])

    out = pd.DataFrame(
        {
            f"{name1}_mean": means1,
            f"{name2}_mean": means2,
            "delta": delta,
            "bf10": bf10,
        },
        index=common,
    )
    return out

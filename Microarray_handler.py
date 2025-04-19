from io import StringIO
from typing import Tuple, Dict, List

import pandas as pd
import pathlib


def read_sdrf(path: str, *,comment_prefix: str = "!",sep: str = "\t",) -> pd.DataFrame:
    path = pathlib.Path(path)
    with path.open("r", encoding="utf-8") as fh:
        content = [line for line in fh if not line.startswith(comment_prefix)]
    df = pd.read_csv(StringIO("".join(content)), sep=sep, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    return df

def _clean_source_name(name: str) -> str:
    """Return the *first* token of a Source Name, trimming whitespace.

    Examples
    --------
    >>> _clean_source_name("GSM765899 1")
    'GSM765899'
    >>> _clean_source_name(" GSM123 ")
    'GSM123'
    """

    return name.strip().split()[0]

def build_methylome_dataframes(
    sdrf_path: str,
    sample_dir: str,
    *,
    group_col: str = "Comment [Sample_source_name]",
    source_col: str = "Source Name",
    value_col: str = "VALUE",
    id_col: str = "Reporter Identifier",
    groups: Tuple[str, str] = ("Newborns", "Nonagenarians"),
    file_suffix: str = "_sample_table.txt",
) -> Tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Return one CpG×samples matrix per group.

    The *Source Name* column is cleaned automatically (``GSM765899 1`` →
    ``GSM765899``) so the corresponding sample table is located.
    """

    sdrf = read_sdrf(sdrf_path)

    mapping = (
        sdrf[[source_col, group_col]]
        .dropna()
        .query(f"`{group_col}` in @groups")
        .assign(clean_id=lambda d: d[source_col].apply(_clean_source_name))
    )
    group_to_samples: Dict[str, List[str]] = (
        mapping.groupby(group_col)["clean_id"].apply(list).to_dict()
    )

    dfs: Dict[str, pd.DataFrame] = {}
    sample_dir = pathlib.Path(sample_dir)

    for group in groups:
        samples = group_to_samples.get(group, [])
        tmp: pd.DataFrame | None = None
        for sample_id in samples:
            file_path = sample_dir / f"{sample_id}{file_suffix}"
            if not file_path.exists():
                print(f"[build_methylome_dataframes] warning: file not found : {file_path}")
                continue

            table = pd.read_csv(file_path, sep="\t", dtype={id_col: str})

            sample_series = (
                table[[id_col, value_col]]
                .rename(columns={value_col: sample_id})
                .set_index(id_col)
                .astype(float)
            )
            tmp = sample_series if tmp is None else tmp.join(sample_series, how="outer")

        dfs[group] = tmp.sort_index() if tmp is not None else None

    return tuple(dfs.get(g) for g in groups)
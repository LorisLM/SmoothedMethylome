![man_savgol17_2](https://github.com/user-attachments/assets/12267599-0e4b-4456-b68c-0f2a28fd22c6)
[![language](https://img.shields.io/badge/language-Python-3776AB)](https://www.python.org/)
[![OS](https://img.shields.io/badge/OS-Linux%20%7C%20Windows%20%7C%20macOS-0078D4)](https://www.python.org/downloads/)
[![CPU](https://img.shields.io/badge/CPU-x86%20%7C%20x64-FF8C00)](https://www.python.org/dev/peps/pep-0011/)
[![dependencies](https://img.shields.io/badge/dependencies-numpy%20%7C%20pandas%20%7C%20scipy%20%7C%20matplotlib%20%7C%20seaborn-brightgreen)](https://pypi.org/)
![release](https://img.shields.io/badge/release_date-April%202025-yellow)

# About

**SmoothedMethylome** is a Python library for denoising methylation signal data through advanced smoothing techniques, including Savitzky-Golay filtering and moving average, followed by peak detection. It is tailored for genomic methylation signal processing, especially in the context of epigenomic studies using CpG methylation profiles.

# Features

- ✅ **Savitzky-Golay smoothing** with customizable window size and polynomial order  
- ✅ **Bootstrap-based noise estimation** to calculate signal-to-noise ratio (SNR)  
- ✅ **Peak detection** from smoothed signals using statistical thresholds  
- ✅ **Visualization utilities** (line plots, violin plots, Manhattan plots)  
- ✅ **Parameter optimization tools** to guide the choice of smoothing parameters  

# Installation

Use the following git command to clone SmoothedMethylome.
```bash
git clone https://github.com/LorisLM/SmoothedMethylome.git
```

# Getting Started

## Data Format

Most functions in the SmoothedMethylome package operate on pandas DataFrames structured around genome-wide methylation statistics. The typical input format includes:

| Column Name               | Type    | Description                                                            |
| ------------------------- | ------- | ---------------------------------------------------------------------- |
| `Chromosome`              | str/int | Chromosome identifier (e.g., `"1"`, `"2"`, `"X"`, `"Y"`)               |
| `Position`                | int     | Genomic coordinate of each CpG probe                                   |
| `Reporter ID` / `CpG ID`  | str     | Unique identifier for the probe (e.g., `cg00000029`)                   |
| `VALUE`                   | float   | Methylation beta value (between 0 and 1) |
| `pvalue` / `log10_pvalue` | float   | Result from a statistical test (e.g., T-test or ANOVA)                 |
| `Group` (optional)        | str     | Sample group or label (e.g., `"Newborns"`, `"Nonagenarians"`)          |

## Signal Smoothing
Raw methylation signals are often noisy due to technical variability or biological dispersion. To enhance interpretability and identify consistent epigenomic patterns, the SmoothedMethylome package includes built-in methods to smooth the signal along the genome, especially across neighboring CpG sites.

This method computes a centered rolling mean over a specified number of CpG positions.
```python
from smoothed_methylome.smoothing import mean_smoothing

smoothed_df = mean_smoothing(df, target_column="log10_pvalue", n=3)
```

This method fits a local polynomial within a sliding window to better preserve peak structure compared to the moving average.
```python
from smoothed_methylome.smoothing import savgol_smoothing

smoothed_df = savgol_smoothing(df, target_column="log10_pvalue", n=11, order=2)
```

> [!NOTE]
> It is recommended to keep the window relatively short (3–7 CpG sites) to reflect biologically realistic co-methylation patterns.

## Use Case
```python
from smoothed_methylome.smoothing import savgol_smoothing
from smoothed_methylome.plotting import plot_manhattan

import pandas as pd

# Data example
df = pd.DataFrame({
    "Chromosome": ["1", "1", "1", "2", "2"],
    "Position": [100000, 200000, 300000, 100000, 200000],
    "log10_pvalue": [0.5, 1.2, 5.3, 0.9, 3.1]
})

df = savgol_smoothing(df, target_column="log10_pvalue", n=3, order=1)

plot_manhattan(df, smooth_column="smooth_result")
```

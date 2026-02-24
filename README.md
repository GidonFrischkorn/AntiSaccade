# Online Supplement: Is the Anti-Saccade Task a Valid Measure of Inhibition?

[![DOI](https://img.shields.io/badge/DOI-10.1037%2Fxge0001808-blue)](https://doi.org/10.1037/xge0001808)

This R package contains the data and code for the online supplement of the paper "Is the Anti-Saccade Task a Valid Measure of Inhibition?" published in the [Journal of Experimental Psychology: General](https://psycnet.apa.org/doi/10.1037/xge0001808).

## Citation

If you use data or code from this package, please cite:

> Frischkorn, G. T. (2025). Is the Anti-Saccade Task a Valid Measure of Inhibition? *Journal of Experimental Psychology: General*. https://doi.org/10.1037/xge0001808

Or in R:

```r
citation("AntiSaccade")
```

## Installation

You can install the R package using the following code:

```r
remotes::install_github("GidonFrischkorn/AntiSaccade")
```

## Analysis Documents

The analysis code is located in the `/docs` folder. All analyses are documented with Quarto documents. You can find the rendered documents [here](https://osf.io/b8hfd/).

**Note:** The Quarto documents load cached Bayesian model fits (`.rds` files) that are not included in this repository due to their size. To re-render the documents, all Bayesian models must be re-fitted, which requires substantial computation time.

## Datasets

If you want to access the documentation of the data sets, type `?dataset_name` in the R console. The data sets included are:

- Experiment 1: `Exp1_data`
- Experiment 2: `Exp2_data`, `Exp2_WMC_data`, and `Exp2_PS_data`
- Supplementary Experiments: `SuppExp1_data`, `SuppExp2_data`, `SuppExp3_data`, `SuppExp4_data`, `SuppExp5_data`, and `SuppExp6_data`
- Reanalysis Data: `Hood_2022_behavioral` and `Unsworth_2023_E2_binomial`

## Functions

The `/R` folder contains two utility functions:

- `ez.dm()` — Calculate EZ-diffusion model parameters from reaction times and accuracy rates. Documentation: `?ez.dm`
- `coeff_alpha()` — Calculate Cronbach's Alpha for a set of measurements. Documentation: `?coeff_alpha`

## SEM Models

The `/SEM` folder contains text files that specify the measurement models estimated for WMC and PS in Experiment 2.

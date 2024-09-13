# Online Supplement: Is the Anti-Saccade Task a Valid Measure of Inhibition?

This R package contains the data and code for the online supplement of the paper "Is the Anti-Saccade Task a Valid Measure of Inhibition?" The paper is currently under review. You can install the R package using the following code:

```r
remotes::install_github("GidonFrischkorn/AntiSaccade")
```

The analysis code is located in the `/docs` folder. All analysis are documented with Quarto documents. You can find the rendered documents [here](https://osf.io/b8hfd/). If you want to access the documentation of the data sets just type `?data_set` in the R console. The data sets included are:

- Experiment 1: `Exp1_data`
- Experiment 2: `Exp2_data`, `Exp2_WMC_data`, and `Exp2_PS_data`
- Supplementary Experiments: `SuppExp1_data`, `SuppExp2_data`, `SuppExp3_data`, `SuppExp4_data`, `SuppExp5_data`, and `SuppExp6_data`

The `/R` folder contains a function to calculate ezDM parameters from reaction times and accuracy rates. Documentation for the function is available by typing `?ez.dm` in the R console.

Finally, the `/SEM` folder contains text files that specify the measurement models estimated for WMC and PS in Experiment 2.

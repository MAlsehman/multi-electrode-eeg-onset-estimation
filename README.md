Estimating EEG Onsets Across Electrodes: A Simulation Study

This repository contains the code, simulated data, and notebooks used in an MSc project investigating the estimation of the onset of a between-condition difference in electroencephalography (EEG) data when multiple electrodes are analysed simultaneously.

The project extends a single-electrode simulation introduced by Guillaume Rousselet (2025) to a nine-electrode setting with fixed spatial weights and a known true effect onset of 160 ms. This makes it possible to measure onset-estimation error directly and to examine how multiple-comparison control, electrode selection, and aggregation across electrodes affect performance.

Simulation source

The construction of the signal and noise in this project is based on the simulation framework used by Rousselet (2025). Components for the signal template, EEG-like noise generation, and power spectrum were used as the basis for data generation, and the simulation was then extended from one electrode to nine electrodes.

The original code and accompanying materials for Rousselet’s simulation are available in his GitHub repository:

https://github.com/GRousselet/onsetsim

Guillaume Rousselet also ran the simulation code used to generate the final data files included in this repository.

Main notebooks

build.notebook.qmd

This notebook explains the construction of the simulation step by step, including the signal template with a known onset at 160 ms, EEG-like noise generation, construction of the two experimental conditions, and application of the fixed spatial weights across the nine electrodes.

sim.notebook.qmd

This is the main analysis notebook. It compares:

* Benjamini-Hochberg (BH)
* Benjamini-Yekutieli (BY)
* Maximum-statistic correction (MAX)
* Change-point detection (CPD)
* MAX-filtered CPD
* Pooled CPD
* Several rules for aggregating electrode-wise onset estimates, including the minimum, quantiles, and median

Performance is evaluated using bias, mean absolute error (MAE), standard deviation (SD), the proportion of estimates occurring before the true onset, and the number of iterations in which no onset was returned.

multiv.notebook.qmd

This notebook provides a multivariate extension of the main analysis. The nine electrodes are combined within Hotelling’s T^2 statistic at each timepoint, and the following strategies are compared:

* Hotelling BH
* Hotelling BY
* Hotelling MAX
* Hotelling CPD

Computational environment

The simulations and statistical analyses were carried out in R.

R scripts

functions.R

Contains the main functions used in the simulations and analyses, including statistical calculations, multiple-comparison procedures, onset-estimation methods, and aggregation of results.

main_sim.R

Runs the main simulation with 50 trials per condition, together with the independent diagnostic simulation used for electrode-level analyses.

size_sim.R

Runs the trial-number sensitivity analysis at 25, 50, 75, and 100 trials per condition.

mv_sim.R

Runs the multivariate simulation used in multiv.notebook.qmd.

analyse_results.R

Contains additional code for organising, checking, and analysing the saved simulation results.

Data files

simres_n50.RData

Results from the main simulation:

10,000 iterations × 50 trials per condition.

simres.RData

Independent diagnostic simulation:

1,000 iterations × 50 trials per condition.

simres_1000_sizes.RData

Results from the trial-number sensitivity analysis:

1,000 iterations at each of 25, 50, 75, and 100 trials per condition.

The 50-trial run contained in this file is independent of the other two 50-trial simulation runs.

simres_n50_mult.RData

Results from the multivariate simulation:

10,000 iterations × 50 trials per condition.

Supporting files in code/

The data-generation workflow uses the following supporting files:

* code/eeg_noise.R for generation of simulated EEG-like noise
* code/erp_template.R for construction of the signal template
* code/meanpower.txt for power-spectrum information used for noise generation
* code/theme_dark_bw.R for the shared plotting theme used across the notebooks

The directory structure should be preserved so that the source() calls and file paths work without modification.

Suggested order for exploring the project

1. build.notebook.qmd to understand how the signal, noise, and multi-electrode data are constructed.
2. sim.notebook.qmd to examine the main comparison of onset-estimation methods.
3. multiv.notebook.qmd to examine the multivariate extension.

To rerun the simulations themselves, use main_sim.R, size_sim.R, and mv_sim.R.

Reference

The simulation framework is based primarily on:

Rousselet, G. A. (2025). Using cluster-based permutation tests to estimate MEG/EEG onsets: How bad is it? European Journal of Neuroscience, 61, e16618.

Full references for the remaining statistical methods are provided in the report and the relevant notebooks.

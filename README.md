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

Standalone R scripts

For ease of access, selected simulation and analysis code from the notebooks is also provided as standalone R scripts. These scripts were created after the main analyses so that the code can be inspected more easily without reading through the full notebooks.

The analyses reported in the project were carried out from the QMD notebooks. The standalone scripts are provided for accessibility and reference and were not used as separate analysis files during the original workflow.

functions.R

Contains extracted helper functions used in the notebook analyses.

main_sim.R

Contains the code corresponding to the main 50-trial simulation and the electrode-level diagnostic simulation.

size_sim.R

Contains the code corresponding to the trial-number sensitivity analysis at 25, 50, 75, and 100 trials per condition.

mv_sim.R

Contains the code corresponding to the multivariate simulation.

analyse_results.R

Contains extracted code for organising, checking, and summarising the saved simulation results.

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

Some simulation components in the code/ directory originate from Guillaume Rousselet’s onsetsim repository and remain subject to the original MIT License. A copy of that licence is included as code/LICENSE_Rousselet.txt.

The directory structure should be preserved so that the source() calls and file paths work without modification.

Suggested order for exploring the project

1. build.notebook.qmd to understand how the signal, noise, and multi-electrode data are constructed.
2. sim.notebook.qmd to examine the main comparison of onset-estimation methods.
3. multiv.notebook.qmd to examine the multivariate extension.

The standalone R scripts can be consulted separately when direct access to the code is preferred over reading the full notebooks.

Reference

The simulation framework is based primarily on:

Rousselet, G. A. (2025). Using cluster-based permutation tests to estimate MEG/EEG onsets: How bad is it? European Journal of Neuroscience, 61, e16618.

Full references for the remaining statistical methods are provided in the report and the relevant notebooks.
# Estimating electroencephalographic onsets across electrodes: a simulation study of how error control and spatial aggregation affect accuracy

Code and data for a simulation study of onset estimation across nine electrodes. The simulation followed the single-electrode framework of Rousselet (2025), retaining the core data-generation structure while extending the design from one to nine electrodes.

MSc Brain Sciences, BIOL5316P, School of Psychology and Neuroscience, University of Glasgow, 2026.

## Quarto notebooks

| Notebook | HTML version | Content |
| -------- | ------------ | ------- |
| `demo.qmd` | [demo](docs/demo.html) | Describes the data-generating model: the temporal signal template, the electroencephalography-like noise, and the spatial weight pattern across the nine electrodes. |
| `simulation.qmd` | [simulation](docs/simulation.html) | Runs the mass-univariate simulation and compares six onset-estimation strategies, with electrode-level diagnostics, aggregation rules, and a trial-number sensitivity analysis. |
| `multivariate.qmd` | [multivariate](docs/multivariate.html) | Runs the multivariate simulation using Hotelling $T^2$ and compares it with the corresponding univariate results. |

Rendered HTML versions of the notebooks are available in the `docs` folder.

## Abbreviations

BH, Benjamini–Hochberg  
BY, Benjamini–Yekutieli  
CPD, change-point detection  
EEG, electroencephalography  
FDR, false discovery rate  
MAE, mean absolute error  
MAX, maximum-statistic procedure  
SD, standard deviation

## Simulation design

Each Monte Carlo iteration represented a complete simulated experiment for one participant, comprising two conditions and nine electrodes. A known true onset of 160 ms was used, with an epoch from 0 to 500 ms sampled every 2 ms. Spatial signal weights were set to `c(0, 0, 0.3, 0.75, 1, 0.75, 0.3, 0, 0)`, so that five electrodes carried the same effect at different amplitudes and four carried no effect. The two conditions compared were noise alone against signal added to noise.

Welch tests were applied at each electrode, together with Benjamini–Hochberg and Benjamini–Yekutieli correction for FDR control, permutation-based maximum-statistic inference, change-point detection using the `changepoint` package, and Hotelling $T^2$ for the multivariate analysis. The random seed was fixed at `set.seed(26)`, and 1,000 permutations were used in the permutation procedures.

No empirical human or animal data were used.

## Dependencies

Additional code dependencies are in the `code` folder. The analyses were carried out in R and the notebooks were written in Quarto, relying on `tidyverse` and `changepoint`. `multivariate.qmd` also uses `DescTools` for an external check of the Hotelling $T^2$ implementation; this check is skipped if the package is not installed.

The notebooks must be run from the root of the repository because the paths to `code` and `data` are relative.

## Simulation results

All simulation results are in the `data` folder, and the figures and reported summaries can be reproduced from the saved results without re-running the simulations.

| File | Simulation |
| ---- | ---------- |
| `simres.RData` | Electrode-level diagnostic, 1,000 iterations, 50 trials per condition |
| `simres_n50.RData` | 10,000 iterations, 50 trials per condition |
| `simres_1000_sizes.RData` | 1,000 iterations at each of 25, 50, 75, and 100 trials per condition |
| `simres_n50_mult.RData` | 10,000 multivariate iterations, 50 trials per condition |

## Code source

Three files from the GitHub repository accompanying Rousselet's study ([GRousselet/onsetsim](https://github.com/GRousselet/onsetsim)) were used without modification: `eeg_noise.R` for electroencephalography-like noise, `erp_template.R` for the temporal signal template, and `meanpower.txt` for the power-spectrum values used to shape the noise. All three are in the `code` folder and retain their original comments indicating this origin.

## Licence

The three files above were released under GPL-3.0. This repository is distributed under the same licence. See `LICENSE`.

## Acknowledgements

This project was supervised by Dr Guillaume Rousselet (University of Glasgow), who also ran the computationally intensive production simulations from the code in this repository.

## Reference

Rousselet, G. A. (2025). Using cluster-based permutation tests to estimate MEG/EEG onsets: How bad is it? *European Journal of Neuroscience*, 61, e16618.

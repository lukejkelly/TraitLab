# TraitLab

## Description

TraitLab implements the Stochastic Dollo model for simulating and fitting binary data on trees using Markov chain Monte Carlo methods.

* Original paper ([Nicholls and Gray, 2008][1])
* Missing data and rate heterogeneity ([Ryder and Nicholls, 2011][2])
* Lateral trait transfer ([Kelly, 2016][3]; [Kelly and Nicholls, 2017][4])
* Couplings to diagnose convergence and construct unbiased estimators ([Kelly, Ryder and Clarté, 2021][5])

The TraitLab [website][6] includes more information, and the [manual][7] included in this repository describes how to conduct analyses. The manual will be updated in September 2021 to include lateral transfer and couplings.

README files in the `borrowing` and `coupling` directories provide more details on these methods.

---

## Requirements

TraitLab runs in Matlab. The lateral transfer and coupling methods use functions in additional Matlab toolboxes; see their README files for details.

Please get in touch if you have any issues.

---

## Analysis

TraitLab reads Nexus-formatted data in a _.nex_ file. Further details are in the [manual][6].

A `startup` file sets up the the path when Matlab is started at the top level of the `TraitLab` directory.

To run an experiment using the GUI:

* Start Matlab in the _TraitLab_ folder
* Execute `TraitLab` to open the analysis GUI

TraitLab will write the settings (.par) and output files to the specified directory.

To analyse samples at the end of a run, open the analysis GUI from the toolbar in the main TraitLab GUI.

Alternatively, experiments can be run in batch mode.
```matlab
batchTraitLab('path to .par file', [optional number to append to output files]);
```
The [manual][6] describes how to run experiments in detail.

---

## Notes

If catastrophes are included in the Stochastic Dollo model and the catastrophe rate `rho` is
* Fixed, then the number of catastrophes a branch of length `t` is `Poisson(rho * t)`
* Allowed to vary, then we integrate `rho` out analytically and the prior number of catastrophes on branches is Negative Multinomial.

MCMC moves on `rho` have been disabled as a consequence; see Chapter 2 of [Kelly (2016)][3] or the supplement of [Kelly and Nicholls (2017)][4] for further details on this calculation

---

## Seeded random numbers

Repeating an experiment with seeded random numbers from the GUI and in batch mode will produce different results due to the order in which components of initialisation are executed. We will fix this issue soon; in the meantime please set the RNG manually before calling TraitLab.

[1]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2007.00648.x/full
[2]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9876.2010.00743.x/full
[3]: https://ora.ox.ac.uk/objects/uuid:6884785c-fccc-4044-b5b2-7a8b7015b2a5
[4]: https://projecteuclid.org/euclid.aoas/1500537738
[5]: https://arxiv.org/pdf/2108.13328.pdf
[6]: http://www.stats.ox.ac.uk/~nicholls/TraitLab/TRAITLAB%20MANUAL.pdf
[7]: https://sites.google.com/site/traitlab/

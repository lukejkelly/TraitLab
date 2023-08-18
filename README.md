# TraitLab

## Description

TraitLab implements the Stochastic Dollo model for simulating and fitting binary data on trees using Markov chain Monte Carlo methods.

* Original paper: [Nicholls and Gray (2008)][1]
* Missing data and rate heterogeneity: [Ryder and Nicholls (2011)][2]
* Lateral trait transfer: [Kelly and Nicholls (2017)][3]
* Couplings to diagnose convergence and construct unbiased estimators: [Kelly, Ryder and Clarté (2021)][4]

See the TraitLab [website][5] for an overview of TraitLab's features and related publications.
The [manual][6] contains a full description of the tools in TraitLab and a step-by-step tutorial.
README files in the `borrowing` and `coupling` directories provide more details on how to run those methods.

---

## Requirements

TraitLab runs in Matlab with its Statistics and Machine Learning toolbox.
Please get in touch if you have any issues.

---

## Analysis

TraitLab reads Nexus-formatted data in a _.nex_ file.

A `startup` file sets up the path when Matlab is started at the top level of the TraitLab directory.

To run an experiment using the GUI:

* Start Matlab in the _TraitLab_ folder.
* Execute `TraitLab` to open the analysis GUI.
* Load data and set the parameters of the experiments then click _start_.

TraitLab will write the settings (.par) and output files to the specified directory.

To analyse samples at the end of a run, open the analysis GUI from the toolbar in the main TraitLab GUI.

Alternatively, experiments can be launched from parameter files.
```matlab
batchTraitLab('<path to .par file>');
```
This allows experiments to be run from the Matlab command window or the terminal.

See the [manual][6] for a full description.


## Example

The _example_ folder contains a synthetic data set and _.par_ files to run single- and coupled-chain experiments.


[1]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2007.00648.x/full
[2]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9876.2010.00743.x/full
[3]: https://projecteuclid.org/euclid.aoas/1500537738
[4]: https://projecteuclid.org/journals/annals-of-applied-statistics/volume-17/issue-2/Lagged-couplings-diagnose-Markov-chain-Monte-Carlo-phylogenetic-inference/10.1214/22-AOAS1676.short
[5]: https://sites.google.com/site/traitlab/
[6]: https://github.com/traitlab-mcmc/TraitLab/blob/master/TRAITLAB_MANUAL.pdf

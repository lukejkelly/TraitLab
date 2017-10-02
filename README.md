# TraitLab with lateral transfer.

## Description

Code which extends the _Stochastic Dollo_ model ([Nicholls and Gray, 2008][1]; [Ryder and Nicholls, 2011][2]) implementation in _[TraitLab][6]_ [(Nicholls, Ryder and Welch)][3] to include lateral trait transfer ([Kelly, 2016][4]; [Kelly and Nicholls, 2017][5]).

---

## Requirements

The code is in _Matlab_ and requires the standard installation with the exception of `bi2de` and `de2bi` from the _Communication Systems_ toolbox. I have included _MEX_ implementations of these functions which can be compiled from within the _TraitLab/_ folder by
```Matlab
addpath Borrowing
mexFiles
```
The compiled functions have the same names and syntaxes as their _Communication Systems_ counterparts.

If you have OpenMP installed then you can compile the functions to run in parallel by uncommenting the relevant parts of _mexFiles.m_, _fastBi2De.c_ and _fastDe2Bi.c_ in _Borrowing/_.

---

## Analysis

To run an analysis:

* Start _Matlab_ and `cd` to the _TraitLab/_ folder.
* Run `TraitLab` to open the analysis GUI
* Set the options for lateral transfer (enable/disable; set initial beta or randomise; allow beta to vary over a MCMC run)
* Proceed according to the instructions in the [manual][3]

**Note:** I will update the code to run in batch mode (where you specify options in a _.par_ file instead) shortly; please let me know if you are interested in this.

To analyse samples at the end of a run, open the analysis GUI from the toolbar in the main _TraitLab_ GUI.

---

## Notes

If catastrophes are included and:
* The option to fix the catastrophe rate, _rho_, at a value is selected, then the number of catastrophes on each branch of the tree _a priori_ is a _Poisson_ random variable
* The option to allow _rho_ to vary is selected, then the number of catastrophes on each branch is a _Negative Binomial_ random variable.

See Chapter 2 of [Kelly (2016)][4] or the supplement of [Kelly and Nicholls (2017)][5] for further details on this calculation

If you opt to start from the true tree in the Nexus file or a tree in an output file, the value for _beta_ in the file (if any) will be used to initialise the MCMC chain, the value in the GUI will be used otherwise.

---

## Synthetic data

The _Nexus_ file for the synthetic data with lateral transfer example (_SIM-B_) in [Kelly (2016)][4] and [Kelly and Nicholls (2017)][5] is included in the _Data/_ folder.


[1]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2007.00648.x/full
[2]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9876.2010.00743.x/full
[3]: http://www.stats.ox.ac.uk/~nicholls/TraitLab/TRAITLAB%20MANUAL.pdf
[4]: https://ora.ox.ac.uk/objects/uuid:6884785c-fccc-4044-b5b2-7a8b7015b2a5
[5]: https://projecteuclid.org/euclid.aoas/1500537738
[6]: https://sites.google.com/site/traitlab/

# Example analyses

We describe how to run coupled and marginal analyses on a synthetic 8-leaf data set

_L8-data.nex_ contains data sampled from the SD model on a random 8-leaf tree.

* root height = 1000, lambda = 0.1, mu = 2.5e-04
* no lateral transfer, missing data or catastrophes

_L8-marginal.par_ is the instructions for a single marginal analysis for `Run_length = 100000` iterations, storing one sample every `Sample_interval = 100` iterations, with no clade constraints (`Impose_clades = 0`) and mu fixed (`Vary_loss_rate = 0`) at its true value (`Initial_loss_rate = 0.221199`, see §7.3 of the manual for details of this parameterisation). We attempt to infer the tree topology (`Vary_topology = 1`) and node times without lateral transfer (`Account_for_lateral_transfer = 0`), catastrophes (`Include_catastrophes = 0`) or missing data  (`Model_missing = 0`).

_L8-coupled.par_ is identical to _L8-marginal.par_ except it runs a pair of chains coupled at lag 20000.
```matlab
Coupled_markov_chains = 1
Coupling_lag = 20000
```
If the meeting time is less than 10^5 iterations then a marginal chain will complete the experiment. If you want the coupled analysis to stop once the chains have met then set `Run_length = 0` in _L8-coupled.par_.


We can run each of these experiments by executing
```matlab
batchTraitLab('example/L8-marginal.par');
batchTraitLab('example/L8-coupled.par');
```
The output will be written to the _example_ folder.

The `Coupling_lag` can also be changed so long as it is an integer multiple of `Sample_interval`. This is a requirement of TraitLab, in that `runmcmcCoupled` calls `MarkovCoupled` which advances the pair of chains by `Sample_interval` iterations before returning the output to `runmcmcCoupled` to write it to file, and so on.

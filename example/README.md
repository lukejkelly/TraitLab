# Example analyses

We describe how to run coupled and marginal analyses on a synthetic 8-leaf data set

_L8-data.nex_ contains data sampled from the SD model on a random 8-leaf tree.

* root height = 1000, lambda = 0.1, mu = 2.5e-04
* no lateral transfer, missing data or catastrophes

_L8-marginal.par_ run a single marginal analysis for 10^5 iterations, storing every 100th sample, with no clade constraints and mu fixed at its true value.  We attempt to infer the tree topology and node times.

_L8-coupled.par_ is identical to _L8-marginal.par_ except it runs a pair of chains at lag 20000. If the coupling time is less than 10^5 iterations then a marginal chain will complete the experiment.

We can run each of these experiments by executing
```matlab
batchTraitLab('example/L8-marginal.par');
```
and
```matlab
batchTraitLab('example/L8-coupled.par');
```
The output will be written to the _example_ folder.

If you want the coupled chain to stop once it has coupled then set `Run_length = 0` in _L8-coupled.par_.

The `Coupling_lag` can also be changed so long as it is an integer multiple of `Sample_interval`. This is a requirement of TraitLab, in that `runmcmcCoupled` calls `MarkovCoupled` which advances the pair of chains by `Sample_interval` iterations before returning the output to `runmcmcCoupled` to write it to file, and so on.

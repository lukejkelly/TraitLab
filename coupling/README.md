# Lagged couplings of phylogenetic MCMC samplers

This code implements the methods described by [Kelly, Ryder and Clarté (2021)][1] to couple pairs of Markov chains targeting a posterior distribution on phylogenetic trees and model parameters. [Jacob et al.][2] describe a generic method for sampling from couplings of MCMC transition kernels, [Biswas et al.][3] use lag-`l` couplings to quantify convergence of Markov chains.

---

## Analysis

The [manual][4] contains a full description of the tools in TraitLab and a step-by-step tutorial for running an analysis.

In order to run a coupled experiment with a pair of chains staggered at lag 100000, edit the following section of the `.par` file from the corresponding marginal analysis.
```
Coupled_markov_chains = 1
% FOLLOWING IS IGNORED WHEN Coupled_markov_chains == 0
Coupling_lag = 100000
```
In order to start the experiment, execute
```matlab
batchTraitLab('path to .par file', [optional digit to add to output file names]);
```
This will create a similar set of output files to the marginal analyses except filenames will have `_x` or `_y` appended to them.
See the README in the `example` directory for an illustration.

As we typically run many pairs of coupled chains, we can use the second argument of `batchTraitLab` to distinguish between them.

Code to generate synthetic data experiments, Slurm submission scripts, and to analyse the resulting output are contained in the [CoupledPhylogenetics](https://github.com/lukejkelly/CoupledPhylogenetics) Github repository.

---

## Algorithm

If `Coupled_markov_chains = 1`, `Coupling_lag = l` and `Run_length = m`, then `batchTraitLab` will call `runmcmcCoupled` which proceeds as follows.
* Create a common initial state.
* Sample `X0` and `Y0` by running short, independent MCMC chains started at the common initial state and which target the prior without catastrophes.
* Advance the `X` chain by `l` iterations of the marginal kernel.
* Do housekeeping on `Y[0]` to match `X[l]` where possible.
* Sample from the coupling of transition kernels followed by housekeeping until `X[t] == Y[t - l]`.
* Write the coupling time `tau` to the `.tau` output file.
* Advance the `X` chain by `max(0, m - tau)` iterations of the marginal kernel.

---

## Code

Many of the functions in `coupling` are coupled versions of the marginal proposal distributions in `core`. Unit tests, some of which are quite slow, compare the marginal distributions of samples from the coupled proposals with their marginal counterparts, and also the proportion of matching samples with what we would expect to occur in theory.

[1]: https://projecteuclid.org/journals/annals-of-applied-statistics/volume-17/issue-2/Lagged-couplings-diagnose-Markov-chain-Monte-Carlo-phylogenetic-inference/10.1214/22-AOAS1676.short?tab=ArticleLinkCited
[2]: https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12336
[3]: https://papers.nips.cc/paper_files/paper/2019/file/aec851e565646f6835e915293381e20a-Paper.pdf
[4]: https://github.com/traitlab-mcmc/TraitLab/blob/master/TRAITLAB_MANUAL.pdf

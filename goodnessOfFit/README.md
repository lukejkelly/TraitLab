# Goodness-of-fit

We can assess model fit via

* Savage–Dickey ratios for time constraints on nodes to assess model fit.
* Bayes factors using the posterior predictive to compare models with and without lateral trait transfer.

## Savage–Dickey ratios

### Relax clade constraint times
[Ryder and Nicholls (2011)][1] consider the support for each internal node constraint.
To relax a constraint, specify it in the `Ignore ages for clades` textbox in the GUI or the corresponding entry in the `.par` file.
See the [manual][2] for an example of relaxing this constraint and also how to generate samples from the corresponding prior in the model without lateral trait transfer.

### Relax leaf constraint times
Each leaf is automatically constrained to be at time 0 unless a clade specifies its sampling range.
[Kelly and Nicholls (2017)][3] consider the support for each leaf node constraint.
To avoid issues with testing on boundaries, they allow the leaf to vary into the past (positive time in TraitLab) and the future (negative time in TraitLab).
For example, in the `L8-data.nex` file in the `example` directory, we could add a clade on taxa 1 to allow its time to vary:
```
CLADE NAME = clade1
ROOTMIN = -10000 ROOTMAX = 1000
TAXA = 1;
```

This approach may be implemented in a future release of TraitLab but for now requires modifying code as described below.

When leaf ages may be negative, we must comment out the following lines in `guifiles/docladeblock.m`:
```matlab
61: if any(clade{n}.rootrange < 0)
62:   % error in reading
63:   badclade = [badclade n];
64: end
```

Allowing leaf times to vary may also vary the tree height calculation.
If we are placing an approximately flat prior on the root age and specify clades that do not have upper bounds, then if we do not want the root time to be affected by an offset leaf we need to edit `core/LogPrior.m` to ignore
```matlab
13: height = s(Root).time - max(min([s.time]), 0);
```

To generate samples from the prior in the SDLT model, the simplest option is to replace
```matlab
[intLogLik, logLkd] = logLkd(state);
```
in `core/logLkd.m` by
```matlab
[intLogLik, logLkd] = deal(0);
```

## Computing Savage–Dickey ratios
To compute Savage–Dickey ratios and plot the output, we use the functions in the _goodnessOfFit_ directory as below.
In this example, suppose we want to perform goodness-of-fit tests by relaxing constraint #3 in the synthetic data example in figures S12 and S13 the supplement of [Kelly and Nicholls (2017)][3].

```matlab
addpath goodnessOfFit;

% Where I kept the output files for my goodness-of-fit analyses
pathToPriorSamples = 'output/SIM_B_P_C3.nex';
pathToSDLTSamples  = 'output/SIM_B_B_C3.nex';
pathToSDSamples    = 'output/SIM_B_N_C3.nex';

% Each output file contains 1 initial state and 30,000 samples (thinned from a
% chain with 3 million iterations); we discard the first 10,001 entries as
% burn-in and analyse the rest
sampInds = 10002:1:30001;

% The node is the most recent common ancestor of leaves 6–10; the age range of
% the constraint is 200 to 500 years before the reference node (at time 0 here,
% see note below for more details)
ancestorOfNodes = {'6', '7', '8', '9', '10'};
constraints     = [200, 500];
referenceNode   = '2';

% Getting sampled node times for each model
p = getConstraintNodeTimes(pathToPriorSamples, ancestorOfNodes, sampInds, referenceNode);
b = getConstraintNodeTimes(pathToSDLTSamples ancestorOfNodes, sampInds, referenceNode);
n = getConstraintNodeTimes(pathToSDSamples, ancestorOfNodes, sampInds, referenceNode);

% Histograms of prior and posterior node times; bin width of 20 years
histogram(b, 'BinWidth', 20, 'Normalization', 'pdf'); hold on
histogram(n, 'BinWidth', 20, 'Normalization', 'pdf');
histogram(p, 'BinWidth', 20, 'Normalization', 'pdf');
plot(repmat(constraints, 2, 1), get(gca, 'YLim'), 'k', 'LineWidth', 2); hold off

alpha(0.5); axis('tight');
legend('SDLT posterior', 'SD posterior', 'SDLT/SD prior');
xlabel('Node time (years before present)', 'Interpreter', 'LaTeX');
ylabel('Relative frequency', 'Interpreter', 'LaTeX')

% Savage–Dickey ratios
propPrior = propInRange(p, constraints);
propSDLT  = propInRange(b, constraints);
propSD    = propInRange(n, constraints);

logSD = log(propPrior) - log([propSDLT, propSD]);

h = stem(1:2, logSD, 'MarkerFaceColor', 'auto'); xlim([0.5, 2.5]);
h.Parent.XTick = 1:2; h.Parent.XTickLabel = {'SDLT', 'SD'};
xlabel('Model', 'Interpreter', 'LaTeX');
ylabel('Log-Savage--Dickey ratio', 'Interpreter', 'LaTeX'); ;

```
We need to specify a reference node for `getConstraintNodeTimes` as when the tree is written to file, the node times times are shifted so that the most recent node has time 0. This is only an issue for unconstrained leaf nodes which can have negative node times as we do not allow ancestral node times to drop below their offspring at time 0.

## Bayes factors using the posterior predictive

*Note:* This may be implemented in a future release of TraitLab but is currently deprecated so the description remains here for reference.

---

To compute Bayes factors using the posterior predictive distribution, we fit the model to one partition of the data and estimate the posterior predictive on the remainder. The simplest way to do this is to construct training and test partitions of the data as separate nexus input files using the `partitionData` function in the `goodnessOfFit` directory. For example:
```matlab
addpath goodnessOfFit;
partitionData('data', 'SIM_B', 2/3);
```
randomly splits `data/SIM_B.nex` into `data/SIM_B-train.nex` and `data/SIM_B-test.nex` with roughly two-thirds of the columns in the full data set going into the training partition.

 Currently, the only way to store the locations of catastrophes along branches in the SDLT model is to save the state of the chain. To do so, use the global variable `SAVESTATES` as a flag:
```matlab
global SAVESTATES; SAVESTATES = 1;
```
When the next MCMC analysis is started, the `saveStates` directory will be created and the `state` struct will be saved every time the summary statistics are written to file. This is not implemented for running couplings of Markov chains.

We use the samples from the posterior on the training data to estimate the posterior predictive distribution for the test data partition.
```matlab
addpath borrowing goodnessOfFit;
lPL = logPredictiveFromStates(pathToData, testData, outFile, sInds, misDat, lostOnes);
```
See `help logPredictiveFromStates` for details of the inputs. `lPL` is an array containing the log predictive likelihood of the test data at each state sampled from the posterior of the model fit to the training data.

To get a more stable estimate of the log marginal predictive likelihood for estimating Bayes factors, use the [log-sum-exp](https://en.wikipedia.org/wiki/LogSumExp#log-sum-exp_trick_for_log-domain_calculations) trick.
```matlab
logMeanX = @(logX) max(logX) + log(mean(exp(logX - max(logX))));
logMeanX(lPL)
```
The log-predictive calculation is based on an unnormalised likelihood. This does not matter when comparing the SDLT and SD models as the functional form of the likelihood is the same in both cases so the unknown constants cancel, it's how the likelihood parameters are calculated that differs.

For example, if we fit the SD and SDLT models to `data/SIM_B-train.nex` created above, with the output file stems set to `SIM_B-train_N` and `SIM_B-train_B` respectively, then we can estimate a Bayes factor as follows.
```matlab
addpath borrowing goodnessOfFit;
lPL_N = logPredictiveFromStates('data', 'SIM_B-test', 'SIM_B-train_N', sInds, misDat, lostOnes);
lPL_B = logPredictiveFromStates('data', 'SIM_B-test', 'SIM_B-train_B', sInds, misDat, lostOnes);
bF = logMeanX(lPL_B) - logMeanX(lPL_N);
```
where `sInds` is chosen by the user, and `misDat` (model missing data) and `lostOnes` (discard singleton patterns) match the settings for the experiments (see `Model_missing` and `Account_rare_traits` in the corresponding .par output files).

[1]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9876.2010.00743.x/full

[2]: https://github.com/traitlab-mcmc/TraitLab/blob/master/TRAITLAB_MANUAL.pdf
[3]: https://projecteuclid.org/euclid.aoas/1500537738
[4]: https://ora.ox.ac.uk/objects/uuid:6884785c-fccc-4044-b5b2-7a8b7015b2a5

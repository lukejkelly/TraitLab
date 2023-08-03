# TraitLab with lateral transfer

We extend the Stochastic Dollo to include lateral trait transfer whereby instances of traits attempt to transfer copies of themselves to other evolving species at rate `beta` ([Kelly, 2016][4]; [Kelly and Nicholls, 2017][5]).

The likelihood parameters are given by the solution of a sequence of initial value problems which grow exponentially in the number of taxa.

---

## Requirements

This method uses the `bi2de` and `de2bi` from the _Communication Systems_ toolbox. If they are not available then _MEX_ implementations of these functions can be compiled from within the _TraitLab_ folder by executing
```Matlab
addpath borrowing 
mexFiles
```
The compiled functions have the same names and syntaxes as their _Communication Systems_ counterparts.

If you have OpenMP then you can compile the functions to run in parallel by uncommenting the relevant parts of _mexFiles.m_, _fastBi2De.c_ and _fastDe2Bi.c_ in _borrowing_.

I've checked the Matlab code on Linux, Mac and PC and the C code on Linux (gcc) and Mac (Clang). Get in touch if you have any issues.

---

## Analysis

To run an analysis with lateral transfer in the GUI:

* Enable/disable lateral transfer
* Initialise `beta` at a fixed value or sample it randomly
* Allow `beta` to vary over a MCMC run

Large values of the parameters slow down the ODE solver to compute the likelihood so we advise against setting high initial values.

Traits recorded in a single taxon are evidence against lateral transfer so uncheck the option 'Account for rare traits' to keep these traits in the analysis.

Alternatively, experiments can be run in batch mode. The full details are in the [manual][3] with the addition of the following lines in the _.par_ file. For example, to account for lateral transfer with a random initial value and explore the corresponding posterior:
```
Account_for_lateral_transfer = 1
% FOLLOWING IS IGNORED WHEN Account_for_lateral_transfer == 0
Vary_borrowing_rate = 1
Random_initial_borrowing_rate = 1
% NEXT LINE IS IGNORED WHEN Random_initial_borrowing_rate == 1
Initial_borrowing_rate = 0.00184681
```

To analyse samples at the end of a run, open the analysis GUI from the toolbar in the main TraitLab GUI. For goodness-of-fit, see the section below.

---

## Notes

When accounting for lateral transfer, we require the catastrophe locations so include their density in the above calculations.

When starting runs from a tree in the output file, say, the initial catastrophe tree will be different as catastrophe locations are not currently stored in the output, only their branch counts. See the below section 'Bayes factors using the posterior predictive' for more details.

If you opt to start from the true tree in the Nexus file or a tree in an output file, the value for _beta_ in the file (if any) will be used to initialise the MCMC chain; the value in the GUI will be used otherwise.

---

## Description of algorithm

With some minor adjustments to the MCMC algorithm (to account for the extra parameter and the fact that catastrophe locations are included in the state), the lateral transfer code in _borrowing_ focuses on the likelihood calculation.

* `borrowingParameters` stores persistent variables for manipulating patterns and calculating expected pattern frequencies and the likelihood.
* `stype2Events` embeds the tree in the plane and reads off the branching, catastrophe and extinction events on the tree in order from the root to the leaves.
* `solveOverTree` calculates the expected pattern frequencies through a sequence of differential equations between events on the tree and updates at branching and catastrophe events.
* `patternCounts` calculates the Multinomial (`lambda` integrated out under improper Gamma prior; Negative Multinomial if Gamma prior is proper) and Poisson likelihoods accounting for systematic removal of site-patterns in the data, if any.

This repository does not implement the Wang--Landau algorithm, parameter approximation scheme or pseudo-marginal method described in [Kelly (2016)][4], but the code is available on request.

---

## Goodness-of-fit

[Ryder and Nicholls, (2011)][2], [Kelly (2016)][4] and [Kelly and Nicholls (2017)][5] perform a number of goodness-of-fit tests:

* Savage–Dickey ratios for node constraints: below.
* Bayes factors using the posterior predictive distribution: below.

### Savage–Dickey ratios
For constraints on **internal** nodes ([Ryder and Nicholls, 2011][2]), it is a case of removing/loosening the time constraint in the input nexus file. For example, for one of the constraints in the synthetic data sets in `data/`, we replace
```
CLADE  NAME = c3
ROOTMIN = 200 ROOTMAX = 500
TAXA = 6 7 8 9 10;
```
by
```
CLADE  NAME = c1
TAXA = 6 7 8 9 10;
```

Each leaf is automatically constrained to be at time 0 unless a clade is defined to specify its sampling range. To avoid issues with testing on boundaries, we allow the leaf to vary into the past (positive time in TraitLab) and the future (negative time in TraitLab). For example, in the Polynesian data experiments in [Kelly and Nicholls (2017)][5], to loosen the constraint on Hawaiian, we added the following clade to the nexus file:
```
CLADE NAME = Hawaiian
ROOTMIN = -10000 ROOTMAX = 1000
TAXA = Hawaiian;
```
If a leaf was already constrained to lie on an interval then we can just loosen the constraint instead.

When nodes ages may be negative, we must comment out the following lines in `guifiles/docladeblock.m`:
```matlab
61: if any(clade{n}.rootrange < 0)
62:   % error in reading
63:   badclade = [badclade n];
64: end
```

Note that allowing leaf times to vary over large ranges can lead to huge variation in the tree height. If we are placing an approximately flat prior on the root age and specify clades that do not have upper bounds, then if we do not want the root time to be affected by an offset leaf we need to edit `core/LogPrior.m` to ignore the corresponding leaf:
```matlab
13: height = s(Root).time - max(min([s.time]), 0);
```
assuming the most recent time of the remaining leaves is 0. Similarly, for a Yule prior then you would need to account for it in calculating `tl` in line 9.

To generate samples from the prior, the simplest option is to replace
```matlab
[intLogLik, logLkd] = logLkd2(state);
```
in `core/logLkd2.m` by
```matlab
[intLogLik, logLkd] = deal(0);
```

To compute Savage–Dickey ratios and plot the output, we use the functions in the _goodnessOfFit_ directory as below. In this example, suppose we want to perform goodness-of-fit tests by relaxing constraint #3 in the synthetic data example above ([Kelly, (2016)][4], bottom-left component of Figure 5.7 and left-hand component of Figure 5.8).

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

### Bayes factors using the posterior predictive
To compute Bayes factors using the posterior predictive distribution, we fit the model to one partition of the data and estimate the posterior predictive on the remainder. The simplest way to do this is to construct training and test partitions of the data as separate nexus input files using goodnessOfFit/partitionData.m. For example:
```matlab
addpath goodnessOfFit;
partitionData('data', 'SIM_B', 2/3);
```
randomly splits `data/SIM_B.nex` into `data/SIM_B-train.nex` and `data/SIM_B-test.nex` with roughly two-thirds of the columns in the full data set going into the training partition.

 Currently, the only way to store the locations of catastrophes is to save the state of the chain. To do so, use the global variable `SAVESTATES` as a flag:
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

---
## Synthetic data

To generate a synthetic data set with lateral transfer, do not use the GUI. Instead, we use the `borrowing/simBorCatDeath` function.

The following example simulates a tree and generates data from the SD with lateral transfer model.

```matlab
global LEAF ROOT

% Generate an isochronous ten-leaf tree with exponential branching rate 0.001
L = 10;
theta = 1e-3;
s = ExpTree(L, theta);

% Adding offset leaves
leaves = find([s.type] == LEAF);
s(leaves(1)).time = rand * s(s(leaves(1)).parent).time;
draw(s);

% Alternatively, one could read in a tree from a nexus file
%   s = nexus2stype('path to nexus file'); draw(s);

% Adding catastrophes
% Catastrophes are identified by their relative locations in `catloc`: a
% catastrophe on with catloc 0.5 occurs at time
%     s(i).time + 0.5 * [s(s(i).parent).time - s(i).time]
% Catastrophes cannot go on the branch linking the root and Adam nodes.
nodes = find([s.type] < ROOT);
s(nodes(1)).cat = 1; s(nodes(1)).catloc = [0.2];
s(nodes(10)).cat = 1; s(nodes(10)).catloc = [0.5];

% Read branching, catastrophe events etc. into struct tEvents and
% corresponding right-left ordering of leaves rl
[tEvents, rl] = stype2Events(s);

% Parameters of SDLT process
tr = [0.1;  ...  % lambda
      2.5e-4; ...  % mu
      2.5e-4; ...  % beta
      1 / 3]; % kappa
xi = fliplr([s(rl).xi]); % Missing data parameters.

% Generate N * L array of binary site patterns, where N is the total number
% of traits generated across the tree and column l of D is entry L + 1 - l
% of rl
D = simBorCatDeath(tEvents, tr);

% Removing empty site-patterns
D = D(sum(D, 2) > 0, :);

% Masking matrix to incorporate missing data
M = (rand(size(D)) > repmat(xi, size(D, 1), 1));
D(M) = 2;

% Adding data to leaves
for l = 1:L
  s(rl(l)).dat = D(:, L + 1 - l)';
end

% Write to Nexus file
sFile = stype2nexus(s, '', 'BOTH', '', '');
fid = fopen(['data', filesep, 'SIM_B_10.nex'], 'w');
fprintf(fid, sFile);
fclose(fid);
```

The code becomes slower as _lambda_, _mu_, _beta_ and the length of the tree increase.


[1]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2007.00648.x/full
[2]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9876.2010.00743.x/full
[3]: http://www.stats.ox.ac.uk/~nicholls/TraitLab/TRAITLAB%20MANUAL.pdf
[4]: https://ora.ox.ac.uk/objects/uuid:6884785c-fccc-4044-b5b2-7a8b7015b2a5
[5]: https://projecteuclid.org/euclid.aoas/1500537738
[6]: https://sites.google.com/site/traitlab/
[7]: https://www.jstor.org/stable/2291091

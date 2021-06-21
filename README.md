# TraitLab with lateral transfer.

## Description

We extend the _Stochastic Dollo_ model ([Nicholls and Gray, 2008][1]; [Ryder and Nicholls, 2011][2]) implementation in [TraitLab][6] [(Nicholls, Ryder and Welch)][3] to include lateral trait transfer ([Kelly, 2016][4]; [Kelly and Nicholls, 2017][5]).

---

## Requirements

The code is in Matlab and requires the standard installation with the exception of `bi2de` and `de2bi` from the _Communication Systems_ toolbox. If `bi2de` and `de2bi` are unavailable then _MEX_ implementations of these functions can be compiled from within the _TraitLab_ folder by
```Matlab
addpath borrowing
mexFiles
```
The compiled functions have the same names and syntaxes as their _Communication Systems_ counterparts.

> If you have OpenMP then you can compile the functions to run in parallel by uncommenting the relevant parts of _mexFiles.m_, _fastBi2De.c_ and _fastDe2Bi.c_ in _borrowing_.

I've checked the Matlab code on Linux, Mac and PC and the C code on Linux (gcc) and Mac (Clang). Get in touch if you have any issues.

---

## Analysis

TraitLab reads Nexus-formatted data in a _.nex_ file. Further details are in the [manual][3].

A `startup.m` file will add the necessary `core`, `guifiles` and `borrowing` directories to the path when Matlab is started at the top level of the `TraitLabSDLT` directory.

To run an experiment using the GUI:

* Start Matlab in the _TraitLabSDLT_ folder
* Run `TraitLab` to open the analysis GUI
* Set the options for lateral transfer (enable/disable; set initial beta or randomise; allow beta to vary over a MCMC run)
* Proceed according to the instructions in the [manual][3]

> Traits recorded in a single taxon are evidence against lateral transfer so uncheck the option _'Account for rare traits'_ to keep these traits in the analysis.

Alternatively, experiments can be run in batch mode. The full details are in the [manual][3] with the addition of the following lines in the _.par_ file. For example, to account for lateral transfer with a random initial value and explore the corresponding posterior:
```
Account_for_lateral_transfer = 1
% FOLLOWING IS IGNORED WHEN Account_for_lateral_transfer == 0
Vary_borrowing_rate = 1
Random_initial_borrowing_rate = 1
% NEXT LINE IS IGNORED WHEN Random_initial_borrowing_rate == 1
Initial_borrowing_rate = 0.00184681
```
Copy the _.par_ file you want to use to the TraitLabSDLT directory root with the filename _batchtlinput.txt_. Alternatively, create a symlink. In the Matlab console, execute `batchTraitLab`, or from Bash do `matlab -r "batchTraitLab; exit"`.

To analyse samples at the end of a run, open the analysis GUI from the toolbar in the main TraitLab GUI.

For goodness-of-fit, see the section below.

In contrast to the [paper][5], time in TraitLab is in _years before the present_ so it _increases into the past_:
* t > 0, the past
* t = 0, the present
* t < 0, the future

This is important when specifying clade constraints. We return to this in the section on goodness-of-fit below.

---

## Notes

If catastrophes are included in the Stochastic Dollo model and the catastrophe rate `rho` is
* Fixed, then the number of catastrophes a branch of length `t` is `Poisson(rho t)`
* Allowed to vary, then we integrate `rho` out analytically and the number of catastrophes on each branch is a Negative Binomial random variable.

When accounting for lateral transfer, we require the catastrophe locations so include their density in the above calculations.

* When starting runs from a tree in the output file, say, the initial catastrophe tree will be different as catastrophe locations are not currently stored in the output, only their branch counts. See the below section 'Bayes factors using the posterior predictive' for more details.

MCMC moves on _rho_ have been disabled as a consequence. See Chapter 2 of [Kelly (2016)][4] or the supplement of [Kelly and Nicholls (2017)][5] for further details on this calculation

If you opt to start from the true tree in the Nexus file or a tree in an output file, the value for _beta_ in the file (if any) will be used to initialise the MCMC chain; the value in the GUI will be used otherwise.

## Seeded random numbers
Repeating an experiment with seeded random numbers from the GUI using batch mode and vice versa will produce different results due to the order in which components of initialisation are executed.

**April 2020:** Some initialisation functions (e.g. `pop` called in `TraitLab`) modify the state of the random number generator before the seed is set in `fullsetup>initMCMC` after the 'Start' button is clicked so until this is fixed just execute `rng(<your seed>);` before executing `TraitLab` or `batchTraitLab`.

---

### Example

The Nexus files for the synthetic data sets in [Kelly (2016, Chapter 4)][4] and [Kelly and Nicholls (2017, supplement)][5] are included in the _data_ folder.

---

## Description of algorithm

With some minor adjustments to the MCMC algorithm (to account for the extra parameter and the fact that catastrophe locations are included in the state), the lateral transfer code in _borrowing_ focuses on the likelihood calculation.

* `borrowingParameters` stores persistent variables for manipulating patterns and calculating expected pattern frequencies and the likelihood.
* `stype2Events` embeds the tree in the plane and reads off the branching, catastrophe and extinction events on the tree in order from the root to the leaves.
* `solveOverTree` calculates the expected pattern frequencies through a sequence of differential equations between events on the tree and updates at branching and catastrophe events.
* `patternCounts` calculates the Multinomial (_lambda_ integrated out) and Poisson likelihoods accounting for systematic removal of site-patterns in the data, if any.

If you are interested in the system for estimating the likelihood parameters in my [thesis][4] and the corresponding pseudo-marginal inference scheme then please get in touch.

---

## Goodness-of-fit

If you are interested in performing the various goodness-of-fit tests described in [Ryder and Nicholls, (2011)][2], [Kelly (2016)][4] and [Kelly and Nicholls (2017)][5]:
* Savage–Dickey ratios for node constraints: below.
* Bayes factors using the posterior predictive distribution: below.
* Posterior exploration via Wang–Landau, etc.: forthcoming, but get in touch in the meantime.

### Savage–Dickey ratios
For constraints on **internal** nodes ([Ryder and Nicholls, 2011][2]), it is a case of removing/loosening the time constraint in the input nexus file. For example, for one of the constraints in the synthetic data sets, I replaced
```
CLADE  NAME = c3
ROOTMIN = 200 ROOTMAX = 500
TAXA = 6 7 8 9 10;
```
by
```
CLADE  NAME = c3
TAXA = 6 7 8 9 10;
```
Alternatively, one could use the option to ignore clade ages in the GUI.

Each **leaf** is automatically constrained to be at time 0 *unless* a clade is defined. Therefore, to relax this automatic constraint, we define a clade which allows the leaf time to vary. To avoid issues with testing on boundaries, we allow the leaf to vary into the past (positive time in TraitLab) and the future (negative time in TraitLab). For example, in the Polynesian data experiments, to loosen the constraint on Hawaiian, I added
```
CLADE NAME = Hawaiian
ROOTMIN = -10000 ROOTMAX = 1000
TAXA = Hawaiian;
```
to the nexus file. If a leaf was already constrained to lie on an interval then we can just loosen the constraint instead.

When nodes times may be negative, one must comment out the following lines in `guifiles/docladeblock.m`:
```matlab
61: if any(clade{n}.rootrange < 0)
62:   % error in reading
63:   badclade = [badclade n];
64: end
```

Note that allowing leaf times to vary over large ranges can lead to huge variation in the tree height. If we are placing an approximately flat prior on the root time and we specify clades that do not have upper bounds (on the node time before present), and we do not want the root time to be affected by an offset leaf, then we need to edit `core/LogPrior.m` to ignore the corresponding leaf:
```matlab
13: height = s(Root).time - max(min([s.time]), 0);
```
assuming the most recent time of the remaining leaves is 0. Similarly, for a Yule prior then you would need to account for it in calculating `tl` in line 9.

To generate samples from the prior, the simplest option is to replace
```matlab
4: [intLogLik, logLkd] = logLkd2_m( state );
```
in `core/logLkd2.m` by
```matlab
4: [intLogLik, logLkd] = deal(0);
```
This is much more straightforward than accounting for it inside the SDLT log-likelihood function. Alternatively, you could set a global variable to switch these options on and off.

To compute Savage–Dickey ratios and plot the output, we use the functions in the _goodnessOfFit_ directory as below. In this example, suppose we want to perform goodness-of-fit by relaxing constraint #3 in the synthetic data example above ([Kelly, (2016)][4], bottom-left component of Figure 5.7 and left-hand component of Figure 5.8).

```matlab
% Set up global variables and workspace
GlobalSwitches; GlobalValues;
addpath core guifiles borrowing goodnessOfFit;

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
p = getConstraintNodeTimes(pathToPriorSamples, ancestorOfNodes, sampInds, ...
                           referenceNode);
b = getConstraintNodeTimes(pathToSDLTSamples,  ancestorOfNodes, sampInds, ...
                           referenceNode);
n = getConstraintNodeTimes(pathToSDSamples,    ancestorOfNodes, sampInds, ...
                           referenceNode);

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
To compute Bayes factors using the posterior predictive distribution, we need to fit the model to one partition and compute the poster predictive for the remaining data. The simplest way to do this is to construct training and test partitions of the data as separate nexus input files using goodnessOfFit/partitionData.m. For example:
```matlab
addpath goodnessOfFit;
partitionData('data', 'SIM_B', 2/3);
```
randomly splits `data/SIM_B.nex` into `data/SIM_B-train.nex` and `data/SIM_B-test.nex` with roughly two-thirds of the columns in the full data set going into the training partition.

 Currently, the only way to store the locations of catastrophes is to save the state of the chain. To do so, use the global variable `SAVESTATES` as a flag:
```matlab
global SAVESTATES; SAVESTATES = 1;
```
When the next MCMC analysis is started, the `saveStates` directory will be created and the `state` struct will be saved every time the summary statistics are written to file.

Having fitted the model using the training partition of the data, we use the samples from the posterior to estimate the posterior predictive distribution for the test data partition.
```matlab
addpath borrowing goodnessOfFit;
lPL = logPredictiveFromStates(pathToData, testData, outFile, sInds, misDat, ...
                              lostOnes);
```
See `help logPredictiveFromStates` for details of the inputs. `lPL` is an array containing the log predictive likelihood of the test data at each state sampled from the posterior of the model fit to the training data.

To get a more stable estimate of the log marginal predictive likelihood for estimating Bayes factors, use the [log-sum-exp](https://en.wikipedia.org/wiki/LogSumExp#log-sum-exp_trick_for_log-domain_calculations) trick.
```matlab
logMeanX = @(logX) max(logX) + log(mean(exp(logX - max(logX))));
logMeanX(lPL)
```
**Note:** The log predictive is based on an unnormalised likelihood; this doesn't matter when comparing the SDLT and SD models as the functional form of the likelihood is the same in both cases, it's how the likelihood parameters are calculated that differs.

For example, if we fit the SD and SDLT models to `data/SIM_B-train.nex` created above, with the output file stems set to `SIM_B-train_N` and `SIM_B-train_B` respectively, then we can estimate a Bayes factor as follows.
```matlab
addpath borrowing goodnessOfFit;
lPL_N = logPredictiveFromStates('data', 'SIM_B-test', 'SIM_B-train_N', ...
                                sInds, misDat, lostOnes);
lPL_B = logPredictiveFromStates('data', 'SIM_B-test', 'SIM_B-train_B', ...
                                sInds, misDat, lostOnes);
bF = logMeanX(lPL_B) - logMeanX(lPL_N);
```
where `sInds` is chosen by the user, and `misDat` (model missing data) and `lostOnes` (discard singleton patterns) match the settings for the experiments (see `Model_missing` and `Account_rare_traits` in the corresponding .par output files).

---
## Synthetic data

To generate a synthetic data set, do not use the GUI, rather use the `simBorCatDeath` function in _borrowing_ as described below.

To generate data from the same process as _SIM-B_ in the _data_ folder, within the main _TraitLab_ folder run:

```matlab
% Set up global variables and workspace
GlobalSwitches; GlobalValues;
addpath core guifiles borrowing;

% Generate an isochronous ten-leaf tree with exponential branching rate 0.001
L = 10;
theta = 1e-3;
s = ExpTree(L, theta);

% Adding offset leaves
leaves = find([s.type] == LEAF);
s(leaves(1)).time = rand * s(s(leaves(1)).parent).time;
draw(s);

% Alternatively, one could read in a tree. For example,
%   s = nexus2stype(['data', filesep, 'SIM_B.nex']); draw(s);

% Adding catastrophes. The catastrophe location gives the relative position of
% a catastrophe a long the branch from s(i).time to s(s(i).parent).time. Do not
% try to put a catastrophe on the branch linking the root and Adam nodes!
nodes = find([s.type] < ROOT);
s(nodes(1)).cat = 1; s(nodes(1)).catloc = [0.2];
s(nodes(10)).cat = 1; s(nodes(10)).catloc = [0.5];

% Read branching, catastrophe events etc. into struct _tEvents_ and
% corresponding right-left ordering of leaves in plane _rl_
[tEvents, rl] = stype2Events(s);

% Parameters of SDLT process
tr = [0.1;  ...  % lambda
      5e-4; ...  % mu
      5e-4; ...  % beta
      0.221199]; % kappa
xi = fliplr([s(rl).xi]); % Missing data parameters.

% Generate _N * L_ array of binary site patterns, where _N_ is the total number
% of traits generated across the tree and the _l_th column of _D_ is the
% _L + 1 - l_th entry of _rl_
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

The code becomes slower as _lambda_, _mu_, _beta_ and the length of the tree increase. Of course the results will be different due to the different missing data probabilities, random number generator seeds, etc.


[1]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2007.00648.x/full
[2]: http://onlinelibrary.wiley.com/doi/10.1111/j.1467-9876.2010.00743.x/full
[3]: http://www.stats.ox.ac.uk/~nicholls/TraitLab/TRAITLAB%20MANUAL.pdf
[4]: https://ora.ox.ac.uk/objects/uuid:6884785c-fccc-4044-b5b2-7a8b7015b2a5
[5]: https://projecteuclid.org/euclid.aoas/1500537738
[6]: https://sites.google.com/site/traitlab/
[7]: https://www.jstor.org/stable/2291091

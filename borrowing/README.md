# Lateral trait transfer

The functions in this directory implement the stochastic Dollo with lateral transfer (SDLT) model whereby instances of traits attempt to transfer copies of themselves to other evolving species at global rate `beta` ([Kelly, 2016][2]; [Kelly and Nicholls, 2017][3]).

---

## Requirements

This method uses the `bi2de` and `de2bi` from the [Communications toolbox](https://uk.mathworks.com/products/communications.html).
If they are not available then _MEX_ implementations of these functions can be compiled from within the _TraitLab_ folder by executing `addpath borrowing; mexFiles` in the command window.
The compiled functions have the same names and syntaxes as their Communications toolbox counterparts.
The functions can be compiled to run in parallel via OpenMP by uncommenting the relevant parts of _mexFiles.m_, _fastBi2De.c_ and _fastDe2Bi.c_.
The C code has been checked on Linux (gcc) and Mac (Clang).
Please get in touch if you have any issues.

---

## Analysis

The [manual][1] describes in detail how to implement the lateral transfer model.

The likelihood parameters the solution of a sequence of initial value problems which grow exponentially in the number of taxa so this method is only feasible for up to 20 taxa.

Large values of the parameters slow down the ODE solver to compute the likelihood parameters so we advise against setting high initial values.

Traits recorded in a single taxon are evidence against lateral transfer so uncheck the option 'Account for rare traits' to keep these traits in the analysis.

See the README in the `goodnessOfFit` directory for techniques to assess goodness-of-fit and compare models.

---

## Notes

When accounting for lateral transfer, we require the catastrophe locations along branches rather than just the number along each branch.
Currently, only the number of catastrophes along each branch are written to output so catastrophe locations on branches are freshly sampled when starting a run from a tree in an output file.

---

## Description of algorithm

The code in this directory focuses on the likelihood calculation.

* `borrowingParameters` stores persistent variables for manipulating patterns and calculating expected pattern frequencies and the likelihood.
* `stype2Events` embeds the tree in the plane and reads off the branching, catastrophe and extinction events on the tree in order from the root to the leaves.
* `solveOverTree` calculates the expected pattern frequencies through a sequence of linear ordinary differential equations between events on the tree and updates at branching and catastrophe events.
* `patternCounts` calculates the log-likelihood.

This repository does not implement the Wang--Landau algorithm, parameter approximation scheme or pseudo-marginal method described in [Kelly (2016)][2], but the code is available on request.

---
## Synthetic data

To generate a synthetic data set from the lateral transfer model, use the `simBorCatDeath` function within the `borrowing` directory.

The following example simulates a tree and generates data from the SDLT model.

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

% Adding some catastrophes
% Catastrophes are identified by their relative locations in `catloc`: a
% catastrophe on branch i halfway between i and its parent occurs at time
%     s(i).time + 0.5 * [s(s(i).parent).time - s(i).time]
% Catastrophes cannot go on the branch leading into the root node
nodes = find([s.type] < ROOT);
s(nodes(1)).cat = 1; s(nodes(1)).catloc = [0.2];
s(nodes(10)).cat = 1; s(nodes(10)).catloc = [0.5];

% Read branching, catastrophe events etc. into struct tEvents and
% corresponding right-left ordering of leaves rl
[tEvents, rl] = stype2Events(s);

% Parameters of SDLT process to generate data
tr = [0.1;  2.5e-4; 2.5e-4; 1 / 3]; % [lambda; mu; beta; kappa]
xi = fliplr([s(rl).xi]); % missing data parameters

% Generate N * L array of binary site patterns, where N is the total number
% of traits generated across the tree and column l of D is the data at leaf
% L + 1 - l of rl
D = simBorCatDeath(tEvents, tr);

% Removing empty site-patterns for traits which did not survive to the leaves
D = D(sum(D, 2) > 0, :);

% Masking missing data
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

The simulation code becomes slower as _lambda_, _mu_, _beta_ and the length of the tree increase.


[1]: https://github.com/traitlab-mcmc/TraitLab/blob/master/TRAITLAB_MANUAL.pdf
[2]: https://ora.ox.ac.uk/objects/uuid:6884785c-fccc-4044-b5b2-7a8b7015b2a5
[3]: https://projecteuclid.org/euclid.aoas/1500537738

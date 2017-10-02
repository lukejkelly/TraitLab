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

I've checked the code on recent versions of _Matlab_ on Linux, Mac and PC. Get in touch if you have any issues.

---

## Analysis

TraitLab reads Nexus-formatted data in a _.nex_ file. Further details are in the [manual][3].

To run an experiment using the GUI:

* Start _Matlab_ and `cd` to the _TraitLab/_ folder.
* Run `TraitLab` to open the analysis GUI
* Set the options for lateral transfer (enable/disable; set initial beta or randomise; allow beta to vary over a MCMC run)
* Proceed according to the instructions in the [manual][3]

**Note:** I will update the code to run in batch mode (where you specify options in a _.par_ file instead) shortly; please let me know if you are interested in this.

To analyse samples at the end of a run, open the analysis GUI from the toolbar in the main _TraitLab_ GUI.

**Goodness-of-fit:** If you are interested in performing the various goodness-of-fit tests described in [Kelly (2016)][4] and [Kelly and Nicholls (2017)][5] using Savage--Dickey ratios, the Wang--Landau algorithm, etc., then do get in touch.

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

To generate a synthetic data set, do not use the GUI, rather use the `SimBorCatDeath` function in _Borrowing_ as described below.

To generate data from the same process as _SIM-B_, within the main _TraitLab_ folder run:

```matlab
% Set up global variables and workspace.
GlobalSwitches; GlobalValues;
addpath guifiles Borrowing;

% Generate an isochronous ten-leaf tree with exponential branching rate 0.001.
L = 10;
theta = 1e-3;
s = ExpTree(L, theta);

% Adding offset leaves.
leaves = find([s.type] == LEAF);
s(leaves(1)).time = rand * s(s(leaves(1)).parent).time;
draw(s);

% Alternatively, one could read in a tree. For example,
%   s = nexus2stype('Data/SIM_B.nex);

% Adding catastrophes. The catastrophe location gives the relative position of
% a catastrophe a long the branch from s(i).time to s(s(i).parent).time.
% Do not try to put a catastrophe on the branch linking the root and Adam nodes!
nodes = find([s.type] < ROOT);
s(nodes(1)).cat = 1; s(nodes(1)).catloc = [0.2];
s(nodes(10)).cat = 1; s(nodes(10)).catloc = [0.5];

% Read branching, catastrophe events etc. into struct _tEvents_ and corresponding
% right-left ordering of leaves in plane _rl_.
[tEvents, rl] = stype2Events(s);

% Parameters of SDLT process
tr = [0.1;  ...  % lambda
      5e-4; ...  % mu
      5e-4; ...  % beta
      0.221199]; % kappa
xi = [s(rl).xi]; % Missing data parameters.

% Generate _N * L_ array of binary site patterns, where _N_ is the total number
% of traits generated across the tree and the _l_th column of _D_ is the
% _L + 1 - l_th entry of _rl_.
D = SimBorCatDeath(tEvents, tr);

% Removing empty site-patterns.
D = D(sum(D, 2) > 0, :);

% Masking matrix to incorporate missing data.
M = (rand(size(D)) > repmat(xi, size(D, 1), 1));
D(M) = 2;

% Adding data to leaves.
for l = 1:L
  s(rl(l)).dat = D(:, L + 1 - l)';
end

% Write to nexus file.
sFile = stype2nexus(s, '', 'BOTH', '', '');
fid = fopen('Data/SIM_B_10.nex', 'w');
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

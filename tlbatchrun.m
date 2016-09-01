function tlbatchrun(startfile)


% MCMC run and sampling parameters
% specfiy total length of MC run
fsu.RUNLENGTH         = RL    ; % integer > 0
% specify intervals at which samples are to be taken
fsu.SUBSAMPLE         = SS    ; % integer > 0 and < runlength

% Random number generator
% to use the , set the random number generator t
fsu.SEEDRAND          = SR    ; % 0 or 1
fsu.SEED              = SE    ; % real number

% File that data is to come from
fsu.DATAFILE          = DFS   ; % 'filename.nex' - include full path if not in current directory

% Specify either Flat or Yule prior for trees
fsu.TREEPRIOR         = TP    ; % 'flat' or 'yule'
% if Flat prior is chosen, set maximum age that root of tree can take
fsu.ROOTMAX           = RM    ; % real number > 0

% Specify initial tree to start from - there are three choices
% randomtree, treefromfile, or truetree
fsu.MCMCINITTREESTYLE = MI    ; %  'rand', 'file', 'true'
fsu.MCMCINITTREEFILE  = MIF   ; % 'filename.nex' - include full path if not in current directory
fsu.MCMCINITTREE      = MT    ; % integer > 0

% Specify inital
fsu.MCMCINITMU        = IM    ;
fsu.MCMCINITTHETA     = IT    ;
fsu.VERBOSE           = VB    ;

%Specify
fsu.OUTFILE           = OF    ;

fsu.LOSTONES          = LO    ;
fsu.LOST              = LT    ;
fsu.LOSSRATE          = LR    ;
fsu.PSURVIVE          = PS    ;


% specify taxa to be excluded from analysis ([] to include all)
fsu.DATAMASK          = []    ;
% specify characters to be excluded from analysis ([] to include all)
fsu.COLUMNMASK        = []    ;

% Specify whether clades shold be imposed or not
fsu.ISCLADE           =  1   ; % 0 or 1
% specify clades to be excluded from analysis ([] to include all)
fsu.CLADE             = []    ;

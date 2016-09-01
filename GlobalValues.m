%% The values in the top part of the file can be modified by the user.

%FIGURE properties
FONTSIZEINFIG=12; %14 for a bigger font

%Missing data properties
XI=0.99;
WRITEXI=1; % Set to 1 to write a file with values of XI.

%HTML properties
IMGSIZE=900; %width (in pixels) of images on HTML pages

% How to treat gaps in a Nexus file
GAP=2; %if GAP==2, gaps will be treated as missing data;  if GAP==0, they will be treated as absence of trait.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CAUTION: DO NOT MAKE ANY CHANGES BELOW THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%MCMCINIT values:
TRUSTART=1;  %MCMC starts from the known true tree (synthetic data only)
EXPSTART=2;  %use ExpTree to generate a random start tree
OLDSTART=3;  %load a start state from a nexus file (the last state)

%run control - edit here to change behaviour
TESTUP=0;    %=1 will check result of every MCMC update for legal state (very slow)
TESTSS=1;    %=1 will check every SUBSAMPLE updates for legal state (harmless)
GATHER=1;    %=1 will save all workspace vars to disk every hour of elapsed CPU time (harmless)

%TODO alter prior tests elsewhere
%TREEPRIOR values:
YULE=1;
FLAT=2;
%Topology prior values:
LABHIST=1; %uniform on labeled histories
TOPO=2; %uniform on topologies

%VERBOSE values
%VERBOSE = 0/1/2 reporting [off/update count/output graphs]
QUIET=0;
COUNT=1;
GRAPH=2;
JUSTS=3;        % draw only the statistics figure
JUSTT=4;        % draw only the tree figure

%MASKING can be on or off
OFF=0;
ON=1;

%when i is a vertex in tree s, allowed s(i).type values:
LEAF=0;
ANST=1;
ROOT=2;
ADAM=3;

%when i is a language and j is a cognate, allowed data.array(i,j) values:
OUT=0; %language i lacks cognate j if s(i).dat(j)==OUT
IN=1;  %language i has cognate j if s(i).dat(j)==IN 
MIX=2; %language i may or may not have cognate j if s(i).dat(j)==MIX


%switches labels to make the MCMC a bit more readable
VARYP=0;  %complain if ~VARYP but move(9)=1
NARROW=1;
WIDE=2;
VARYNU=0;
DEPNU=1;
VARYLAMBDA=0;

%OTHER: for tree s, used with s.sibling to get "other" child
%if s(i).child=[j,k] then s(j).sibling=1 and s(k).sibling=2
%so we can write k=s(i).child(OTHER(s(j).sibling))
OTHER=[2,1];

%TODO: lost globals to check
%za zb zc C N NS L LEAVES NODES  DEL DELDEL cmove  sTR RootTR LEAVEStr NODEStr ALPHA BETA 
%THETA LOST LOSTONES LANGOUT MAXHEIGHT

%DATASOURCE values: (this variable was previously called 'ACTION')
NEXUS=1;
BUILD=2;

%SYNTHSTYLE values
OLDTRE=1;
NEWTRE=2;

% draw() node text values:
%LEAF
%ANST
COGS=LEAF+ANST+1; %ie COGS is a whole number not equal LEAF or ANST
CONC=COGS+1;

% GUI related variables
STOPRUN = 0; % used as flag to stop MCMC run. Value changed by GUI



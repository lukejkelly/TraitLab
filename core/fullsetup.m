function [data,model,state,output,mcmc] = fullsetup(fsu)

%header file defines globals

global BUILD TESTSS
%setup the mcmc control variables
TWOSEQS = (~fsu.MASKING & fsu.NUMSEQ==2) | (fsu.MASKING & fsu.NUMSEQ-length(fsu.DATAMASK)==2);
mcmc=initMCMC(fsu,TWOSEQS);
%initialise observation and prior models
model=initMODEL(fsu);
%setup the data source control variables
[data.initial,data.true]=initDATA(fsu);
%build or load the data
[data.content,data.true]=makeDATA(model.prior,data.initial,data.true,fsu.GUICONTENT);

if nargout<=2, return; end
%TODO this is not good RJR 25/07/07
if data.initial.source==BUILD && data.initial.synthclades
    model.prior.clade=synthclades(data.true.state.tree,data.initial.numclades,data.initial.originrootbooth,data.initial.accuracy);
end
%generate the start-state for the MCMC
[state,data.true.state]=initSTATE(mcmc,model,data,fsu.MCMCINITTREE);
%is the start state a legal state?
if TESTSS && check(state,data.true.state)
   disp('initial state doesnt check');
   keyboard;pause;
end
%write the start-state to the sample bag
output=initOUTPUT(state,mcmc.update,fsu);

% Luke 01/04/21: adjust mcmc.update (set in initMCMC) if all leaf times are
% fixed (only set afterwards in initState)
if all([state.tree.leaf_has_timerange] == 0)
    mcmc.update.move(11) = 0;
    mcmc.update.move = mcmc.update.move ./ sum(mcmc.update.move);
    mcmc.update.cmove = cumsum(mcmc.update.move);
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mcmc=initMCMC(fsu,TWOSEQS)

global GATHER
global VARYRHO VARYMU VARYLAMBDA VARYKAPPA VARYNU DEPNU VARYP MISDAT MCMCCAT DONTMOVECATS BORROWING VARYBETA GF_START GF_APPROX

VARYRHO=fsu.VARYRHO*fsu.MCMCCAT;
VARYMU=fsu.VARYMU;
VARYKAPPA=fsu.VARYKAPPA*fsu.MCMCCAT; % LUKE 05/10/2013
MCMCCAT=fsu.MCMCCAT;
DONTMOVECATS=fsu.DONTMOVECATS;
MISDAT=fsu.MISDAT;
BORROWING = fsu.BORROWING;
VARYBETA = fsu.VARYBETA && fsu.BORROWING;

% Green's functions things - Luke 12/12/2014
GF_START = fsu.GF_START;
GF_APPROX = fsu.GF_APPROX;

if DEPNU
    fsu.MCMCINITNU=fsu.MCMCINITKAPPA*fsu.MCMCINITLAMBDA/fsu.MCMCINITMU;
else
    disp('Warning: DEPNU is set to 0 (possibly in GlobalValues.m), but TraitLab does not currently support this.');
end

%START STATE
mcmc.monitor=struct('on',{0},'filename',{'dolloProfile'},'data',{[]});
mcmc.initial=struct('seedrand',{fsu.SEEDRAND},'seed',{fsu.SEED},'setup',{fsu.MCMCINITTREESTYLE},'oldstatefile',{fsu.MCMCINITTREEFILE},'mu',{fsu.MCMCINITMU},'p',{fsu.MCMCINITP},'theta',{fsu.MCMCINITTHETA},'nu',{fsu.MCMCINITNU},'rho',{fsu.MCMCINITRHO},'kappa',{fsu.MCMCINITKAPPA},'lambda',{fsu.MCMCINITLAMBDA},'missing',{fsu.MCMCMISS}, 'beta', {fsu.MCMCINITBETA});

if ~MCMCCAT
    mcmc.initial.rho=0;
    mcmc.initial.kappa=0;
    mcmc.initial.nu=0;
end

if ~BORROWING, mcmc.initial.beta = 0; end

%MOVE TYPES
%1 slide
%2 exchange local
%3 exchange wide
%4 balding local
%5 balding wide
%6 random scale
%7 scale subtree
%8 random walk mu (log scale)
%9 random walk p (random power)
%10 missing data
%11 vary leaves
%12 rescale tree above clade bounds
%13 add a catastrophe
%14 delete a catastrophe
%15 random walk rho (log scale)
%16 random walk kappa (log scale)
%17 random walk lambda (log scale)
%18 Move catastrophe to neighbour
%19 Vary XI for one leaf
%20 Vary XI for all leaves
%21 random walk beta (log scale) % Luke
if TWOSEQS      % Luke 23/04/2014 added Beta. % Luke 02/10/2017 removed move on Rho.
   move=[1, 0, 0, 0, 0, 1, 0, VARYMU, VARYP, 0, 0, 1*fsu.ISCLADE, MCMCCAT, MCMCCAT, 0, VARYKAPPA, VARYLAMBDA, fsu.MCMCCAT, fsu.MCMCMISS, fsu.MCMCMISS, VARYBETA]; % Luke 24/10/14 changed move 10 from MISDAT to 0.
   %     1  2  3  4  5  6  7  8       9      10 11 12             13       14       15 16         17          18           19            20            21
   %move=[0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1];
else
   move=[1, 3*fsu.MCMCVARYTOP, fsu.MCMCVARYTOP, 3*fsu.MCMCVARYTOP, fsu.MCMCVARYTOP, 3, 1, 3*VARYMU, VARYP, 0, 3, 3*fsu.ISCLADE, fsu.MCMCCAT, fsu.MCMCCAT, 0, fsu.VARYKAPPA, VARYLAMBDA, fsu.MCMCCAT, fsu.MCMCMISS, fsu.MCMCMISS, VARYBETA]; % Luke - added #21 VARYBETA.
   %     1  2                  3                4                  5                6  7  8         9      10 11 12             13           14           15 16             17          18           19            20            21
   if fsu.DONTMOVECATS %for punctuational bursts.
       move(13)=0;
       move(14)=0;
       move(18)=0;
   end

   %move=[1 0 0 0 0 0 1 VARYMU VARYP MISDAT 1];
   %move=[1 0 0 0 0 1 0 0 0 0 0 0 1 1 1 0 0 0];
   %move=[zeros(1,12),1,1,zeros(1,3),1];
end
move=move./sum(move);
del=0.5;

mcmc.update=struct('Nmvs',{length(move)},'move',{move},'cmove',{cumsum(move)},'del',{del},'deldel',{(1/del)-del},'theta',{fsu.MCMCINITTHETA});
mcmc.runlength=fsu.RUNLENGTH;
mcmc.subsample=fsu.SUBSAMPLE;
mcmc.gather=GATHER;

%% Luke 3/4/20: Replacing old RNG format
%% TODO: Move RNG seeding to start of code
if mcmc.initial.seedrand
    warning(['Setting specified RNG seed is too far into the code as some', ...
             'of the initialisation functions (e.g. pop) called before', ...
             'this use the RNG']);
    % use given seed
    rng(fsu.SEED);
else
    % reset generator to random state
   rng('shuffle'); % rng(sum(100 * clock));
end

if mcmc.monitor.on
   profile on;
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model=initMODEL(fsu)

global YULE FLAT

global LOSTONES; %TODO FIX THIS
LOSTONES=fsu.LOSTONES;

model=pop('model');

switch fsu.TREEPRIOR
case YULE
   model.prior.type=YULE;
   model.prior.rootmax=[];
case FLAT
   model.prior.type=FLAT;
   model.prior.rootmax=fsu.ROOTMAX;
end

model.prior.topologyprior=fsu.TOPOLOGYPRIOR;

model.observe.LostOnes=fsu.LOSTONES;
model.observe.lossrate=fsu.LOSSRATE;
model.observe.p=fsu.PSURVIVE;

model.prior.cat=fsu.MCMCCAT;

model.prior.isclade=fsu.ISCLADE;
if model.prior.isclade
    model.prior.clade=fsu.CLADE;
    model.prior.strongclades=fsu.STRONGCLADES;
    n=size(model.prior.clade,2);
    model.prior.isupboundadamclade=logical([]);
    model.prior.upboundclade=[];
    for k=1:n
        if (~isempty(model.prior.clade{k}.adamrange) && ~isempty(model.prior.clade{k}.rootrange))
            if (model.prior.clade{k}.adamrange(1)<model.prior.clade{k}.rootrange(1))
                model.prior.clade{k}.adamrange(1)=model.prior.clade{k}.rootrange(1);
                warning(['Setting the lower bound on origin in clade ',model.prior.clade{k}.name,' equal to the higher lower bound on root.']);
            end
            if (model.prior.clade{k}.rootrange(2)>model.prior.clade{k}.adamrange(2))
                model.prior.clade{k}.rootrange(2)=model.prior.clade{k}.adamrange(2);
                warning(['Setting the upper bound on root in clade ',model.prior.clade{k}.name,' equal to the lower upper bound on origin.']);
            end
        end
        if (~isempty(model.prior.clade{k}.adamrange) && ((model.prior.type==FLAT && model.prior.clade{k}.adamrange(2)<model.prior.rootmax) || (model.prior.type==YULE && ~isinf(model.prior.clade{k}.adamrange(2)))))
            model.prior.upboundclade=[model.prior.upboundclade,k];
            model.prior.isupboundadamclade=[model.prior.isupboundadamclade,true];
            model.prior.clade{k}.lowlim=model.prior.clade{k}.adamrange(1);
        else
            if (~isempty(model.prior.clade{k}.rootrange) && ((model.prior.type==FLAT && model.prior.clade{k}.rootrange(2)<model.prior.rootmax)) || (model.prior.type==YULE && ~isinf(model.prior.clade{k}.rootrange(2))))
                model.prior.upboundclade=[model.prior.upboundclade,k];
                model.prior.isupboundadamclade=[model.prior.isupboundadamclade,false];
                model.prior.clade{k}.lowlim=model.prior.clade{k}.rootrange(1);
            end
        end
    end
else
    model.prior.clade={};
    model.prior.strongclades=[];
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [initial,true]=initDATA(fsu)

global NEXUS BUILD OLDTRE NEWTRE OFF

initial=pop('initial');

initial.source=fsu.DATASOURCE;
initial.file=fsu.DATAFILE;
initial.synth=fsu.SYNTHSTYLE;
initial.treefile=fsu.SYNTHTREFILE;
initial.tree=fsu.SYNTHTRE;
initial.knowcats=fsu.KNOWCATS;
initial.masking=fsu.MASKING;
initial.mask=fsu.DATAMASK;
initial.is_column_mask=fsu.ISCOLMASK;
initial.column_mask=fsu.COLUMNMASK;
initial.lost=fsu.LOST;


switch initial.source
    case NEXUS
        if fsu.DATASYN == OFF
            true=[];
        else
            % data in nexus file has been synthesized - have true state
            true = fsu.GUITRUE;
        end
    case BUILD
        initial.borrow=fsu.BORROW;
        initial.localborrow = fsu.LOCALBORROW;
        initial.maxdist = fsu.MAXDIST;
        initial.polymorph=fsu.POLYMORPH;
        initial.nmeaningclass=fsu.NMEANINGCLASS;
        initial.rhbranchvar = fsu.LOSSRATEBRANCHVAR;
        initial.rhclassvar = fsu.LOSSRATECLASSVAR;
        initial.missing=fsu.SYNTHMISS;
        initial.synthclades=fsu.SYNTHCLADES;
        initial.numclades=fsu.NUMSYNTHCLADES;
        initial.accuracy=fsu.SYNTHCLADESACCURACY;
        initial.originrootbooth=fsu.SYNTHCLADESTIMES;
        initial.cats=fsu.MCMCCAT;

        true=pop('true');
        true.mu=DeathRate(fsu.LOSSRATE);
        true.p=fsu.PSURVIVE;
        true.br=fsu.BORROWFRAC*true.mu;
        true.beta = fsu.BETA; % LUKE 05/10/2013
        true.vocabsize=fsu.VOCABSIZE;
        true.lambda=true.vocabsize*true.mu;
        true.theta=fsu.THETA;
        true.rho=fsu.RHO;
        true.kappa=fsu.KAPPA;
        if fsu.DEPNU
            true.nu=true.lambda*true.kappa/true.mu;
        else
            true.nu=fsu.NU;
        end


        switch initial.synth
            case OLDTRE
                true.NS=[];
            case {NEWTRE}
                true.NS=fsu.NUMSEQ;
        end


        if ~fsu.MCMCCAT
            true.kappa=0;
            true.rho=0;
            if initial.synth==NEWTRE, true.cat=zeros(1,2*true.NS); end
        end
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [content,true]=makeDATA(prior,initial,true,GUICONTENT)
%global dataold
global NEXUS BUILD

switch initial.source
case NEXUS
    disp(sprintf('\n***Loading data presented from GUI'));
    content = GUICONTENT;
    %content = ObserveData(content,[],[],0);
    if ~isempty(true)
       true.state=makestate(prior,true.mu,true.lambda,true.p,true.rho,true.kappa,content,true.state.tree, [], true.beta);
    else
       true.state=[];
    end
case BUILD
    [true,content]=truthN(prior,initial,true);
%    true=dataold.true;
%   content=dataold.content;
    % make sure we have the specified value of lambda in true.state
    true.state.lambda = true.lambda; %DW 9/1/2007

end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output=initOUTPUT(state,update,fsu)
global ANST LEAF

%OUTPUT STATISTICS
output = pop('output');

%VERBOSE options {QUIET, COUNT, GRAPH} runtime output verbosity
output.verbose = fsu.VERBOSE;

% fsu.OUTFILE and fsu.OUTPATH are the file and path where output are stored.
% Note that the tree output is saved to OUTFILE.nex while stats are saved to OUTFILE.txt
output.file = fsu.OUTFILE;
output.path = fsu.OUTPATH;

%Write the start state of the Markov chain to output
stats = [ state.logprior; state.loglkd; state.tree(state.root).time; state.mu; state.p; 0; state.lambda; state.kappa ; state.rho ; state.ncat ; state.fullloglkd; state.beta]; % Luke 19/01/2014 added 'state.beta' to end.
output.stats = stats;

%proportion accepted - pad with zeros -> size(output.pa,2)=size(output.stats,2)
output.pa = zeros(update.Nmvs,1);

for i=1:2*state.NS
    if any(state.tree(i).type==[ANST LEAF]), state.tree(i).cat=rand(state.cat(i),1)*(state.tree(state.tree(i).parent).time-state.tree(i).time)+state.tree(i).time; end
end

output.trees{1}=wnextree(state.tree,state.root);

%keep track of the number of samples in the bag: output.Nsamp=size(output.stats,2)
output.Nsamp=1;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

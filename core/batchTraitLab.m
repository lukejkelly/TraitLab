function batchTraitLab(run_file, output_file_name_app)
% Run marginal MCMC sampler for input run_file and (if specified) append
% output_file_name_app to output filenames

GlobalSwitches
GlobalValues
addpath('core') % Luke 05/10/2017
addpath('guifiles') %commented out GKN Feb 08; added back in RJR�16�Mar 2011

rng('shuffle');

% Clear persistent variables in SDLT code. LUKE 24/3/20
clear logLkd2_m patternCounts patternMeans

a = readrunfile(run_file);
for i = 1:length(a{1})
    switch a{1}{i}
        case {'Data_file_name','Output_file_name','Output_path_name'}
            % strings which must have value
            if isempty(a{2}{i})
                error('Must have non-empty value for %s',a{1}{i})
            end
            eval([a{1}{i} '= ''' a{2}{i} ''';']);
        case {'Tree_file_name','Omit_clade_list','Omit_clade_ages_list'}
            % strings that could be empty
            eval([a{1}{i} '= ''' a{2}{i} ''';']);
        case {'Use_tree','Omit_taxa_list','Omit_trait_list','Max_root_age','With_seed'}
            % possibly empty values
            if isempty(a{2}{i})
                eval([a{1}{i} '= [];']);
            else
                eval([a{1}{i} '=[' a{2}{i} '];']);
            end
        otherwise % all compulsory numeric values
            if isempty(a{2}{i})
                error('Must have non-empty value for %s',a{1}{i})
            end
            eval([a{1}{i} '=[' a{2}{i} '];']);
    end
end

% check that output file has no appended extension
x = strfind(Output_file_name,'.nex');
if ~isempty(x) && x == length(Output_file_name)-3
    disp('Discarding .nex extension from Output_file_name')
    Output_file_name = Output_file_name(1:x-1);
end

% Add appendix to output file name when doing multiple runs
if nargin == 2
    Output_file_name = [Output_file_name, '-', num2str(output_file_name_app)];
end

% make sure dont overwrite old run
if exist([Output_path_name Output_file_name '.nex'],'file') && ~strcmp('tloutput',Output_file_name)
    ok = input(['\n' Output_file_name ' already exists in directory ' Output_path_name '\nType 1 to continue and overwrite it\n> '],'s');
    if ~strcmp(ok,'1')
        error('Run aborted to avoid deleting old run')
        disp('Above is untrue if run from the command line and you can read this')
    end
end

% try reading the data file
[truetree,content,data.true,data.clade] = nexus2stype(Data_file_name);
if isempty(content.array)
    error('No data found in data file %s',Data_file_name)
else
   L=size(content.language,1);
   disp(sprintf('%5s %12s%% %s', 'index', 'missing data', 'language'))
   for k=1:L,
       disp(sprintf('%5g %12d%% %s', k, round(100*mean(content.array(k,:) == MIX)), content.language{k}))
   end
end
GC = content;
DSN = ~isempty(truetree);
[NS,L] = size(content.array);
if DSN
    GT = data.true;
    GT.state.tree = truetree;
    GT.NS = length(GT.state.tree)/2;
else
    GT = [];
end

% see which prior we have
if ~exist('Yule_prior_on_tree','var')
    Yule_prior_on_tree = 0;
end
if ~exist('Flat_prior_on_tree','var')
    Flat_prior_on_tree = 0;
end
prior = [Yule_prior_on_tree,Flat_prior_on_tree];
if all(prior == 1) || all(prior ==0 )
    error('Exactly one of Flat_prior_on_tree and Yule_prior_on_tree must be 1')
end
if Yule_prior_on_tree == 1
    TP = YULE;
    RM = 0;
else
    TP = FLAT;
    RM = Max_root_age;
    if RM <= 0
        error('Need postive value for  Max_root_age')
    end
end

%GKN 18 Mar 2011 TOPO prior - added topo option to batch
% see which tree topology prior we have
if ~exist('Uniform_prior_on_tree_topologies','var')
    Uniform_prior_on_tree_topologies = 0;
end
if ~exist('Uniform_prior_on_labelled_histories','var')
    Uniform_prior_on_labelled_histories = 0;
end
topoprior = [Uniform_prior_on_tree_topologies,Uniform_prior_on_labelled_histories];
if all(topoprior == 1) || all(topoprior ==0 )
    error('Exactly one of Flat_prior_on_tree and Yule_prior_on_tree must be 1')
end

    if Uniform_prior_on_tree_topologies==1
        TOPOLOGYPRIOR=TOPO;
    else
        TOPOLOGYPRIOR=LABHIST;
    end

VARYMU = Vary_loss_rate;
if Random_initial_loss_rate
    LR = rand;
else
    LR = Initial_loss_rate;
end
IM = DeathRate(LR);

% Lateral transfer LJK 23/3/20
BORROWING = Account_for_lateral_transfer;
if BORROWING
    VARYBETA = Vary_borrowing_rate;
    ISBETARANDOM = Random_initial_borrowing_rate;
    if ISBETARANDOM
        % Same initialisation as startbutt_Callback.m
        MCMCINITBETA = IM * (0.5 + rand);
    else
        MCMCINITBETA = Initial_borrowing_rate;
    end
else
    VARYBETA = 0;
    MCMCINITBETA = 0;
    ISBETARANDOM = 0;
end

 % Parameters for catastrophes
if ~exist('Include_catastrophes','var')
   MCT = 0;
else
    MCT=Include_catastrophes;
end

% Initialize catastrophe parameters
RMIK=0;
MIK=0;
RMIR=0;
MIR=0;

if MCT
    if exist('Random_initial_cat_death_prob','var') && Random_initial_cat_death_prob
        MIK= 0.25 + 0.75 * rand;
        RMIK=1;
    elseif exist('Initial_cat_death_prob','var')
        MIK=Initial_cat_death_prob;
    else
        warning('Catastrophes included, but no value set for kappa (probability of death at a catastrophe). Setting at 0.2.')
        MIK=.2;
    end

    if exist('Random_initial_cat_rate','var') && Random_initial_cat_rate
        MIR=rand;
        RMIR=1;
    elseif exist('Initial_cat_rate','var')
        MIR=Initial_cat_rate;
    else
        warning('Catastrophes included, but no value set for rho (rate of catastrophes). Setting at 0.002.')
        MIR=.002;
    end
end

if exist('Model_missing','var')
    MISS = Model_missing;
else
    warning('Parameter Model_missing is not defined. Setting to 1.')
    MISS=1;
end

% find what kind of initial tree we have
% define MI, IT, MIF, MT, TN
TN = 0;
IT =0;
MIF = '';
MT = [];

if ~exist('Start_from_rand_tree','var')
   Start_from_rand_tree = 0;
end
 if ~exist('Start_from_tree_in_output','var')
   Start_from_tree_in_output = 0;
end
if ~exist('Start_from_true_tree','var')
   Start_from_true_tree = 0;
end

vals = [Start_from_rand_tree,Start_from_tree_in_output,Start_from_true_tree];
if sum(vals) ~= 1
    error('Exactly one of Start_from_rand_tree, Start_from_tree_in_output, Start_from_true_tree must be 1')
end

switch find(vals)
    case 1
        % use Exptree to generate random initial tree
        MI = EXPSTART;
        IT = Theta;
    case 2
        % use a tree stored in a nexus output file to start
        MI = OLDSTART;
        MIF = Tree_file_name;
        if isempty(MIF)
            error('Start_from_tree_in_output but Tree_file_name not specified')
        else
            TN = Use_tree;
            % load tree from output file - leave off extension when calling
            [oldoutput,ok]=readoutput(MIF(1:end-4));
            if ok && TN <= length(oldoutput.trees)
                % make sure that there
                MT = pop('state');
                MT.tree = rnextree(oldoutput.trees{TN});
                warning('Parameters values from same output as tree');
                %GKN 18/3/11 alter whatever we concluded above
                MT.mu = oldoutput.stats(4,TN);
                MIK = oldoutput.stats(8,TN);
                MIR = oldoutput.stats(9,TN);
                IM=MT.mu;
                MT.rho=MIR;
                MT.p=1;
                MT.kappa=MIK;

                %Add catastrophes %GKN 18/3/11 was '.cat.nex'
                if MCT && exist([MIF(1:end-4) 'cat.nex'],'file')
                    scat=rnextree(oldoutput.cattrees{TN});
                    [MT.tree MT.cat]=CatTreeToList(scat,MT.tree);
                else
                    MT.cat = [];
                end
                % estimate theta at 1/E[edge length]
                IT = (length(MT.tree)-2)/TreeLength(MT.tree,find([MT.tree.type]==ROOT));

                % Lateral transfer LJK 23/3/20
                if BORROWING
                    % We initialise beta = 0 in pop('state').
                    MT.beta = oldoutput.stats(12, TN);
                    warning('Catastrophe locations different as not stored in .cat file');
                end
            else
                error('\nProblem loading initial tree from file.\nCheck Tree_file_name (%s) and Use_tree number (%g).\nNote that file name must include .nex extension.',MIF,TN)
            end
        end
    case 3
        % use true tree to start from
        MI = TRUSTART;
        IT = Theta;
        MT = pop('state');
        MT.tree = truetree;
        MT.mu = data.true.mu;
        MT.p = data.true.p;
        if IT == 0
            % estimate theta at 1/E[edge length]
            IT = (length(MT.tree)-2)/TreeLength(MT.tree,find([MT.tree.type]==ROOT));
        end
        if ~isempty(MT.mu) && LossRate(data.true.mu)>0 && Vary_loss_rate
            % specified a starting value to vary mu from
            IM=MT.mu;
            disp(sprintf('Ignoring the fixed trait death rate set as starting value (%g)',LossRate(data.true.mu)));
            disp(sprintf('Using the trait death rate %g imported with the true tree to initialise MCMC',LossRate(data.true.mu)));
        end
        MT.beta = data.true.beta;
        if BORROWING && MT.beta > 0
            fprintf('Ignoring the trait transfer rate set as starting value (%g) and using\n', MCMCINITBETA);
            fprintf('the rate %g imported with the true tree to initialise the MCMC instead', MT.beta);
            MCMCINITBETA = MT.beta;
            warning('Catastrophe locations different as not stored in .cat file')
        end
end

IC = Impose_clades;
if IC
    CL = data.clade;
    if ~isempty(data.clade)
        disp(sprintf('%1.0f clades found in file',length(data.clade)));
        for hdci=1:length(data.clade)
            disp([sprintf('Clade %3d ',hdci),sprintf(' %s',data.clade{hdci}.language{:})]);
        end
    end
else
    CL = [];
end
if ~exist('Omit_clade_ages_list'), Omit_clade_ages_list=''; end  % For backwards compatibility RJR 07/06/11
[ok,IC,CL,CLM,CLAM]=initclades(Omit_clade_list, IC, CL, Omit_clade_ages_list);

 % these parameters are arbirtary
 ST = NEWTRE;
 STF = '';
 STR = [];
 VS = 195;
 TH = 1/1000;
 BW = ON;
 BF = 0.15;
 PS=1;
 IP=1;




%write the control variables into structures used by fullsetup
fsu=pop('fsu');
fsu.RUNLENGTH         = Run_length    ;
fsu.SUBSAMPLE         = Sample_interval    ;
fsu.SEEDRAND          = Seed_random_numbers ;
fsu.SEED              = With_seed    ;
fsu.DATASOURCE        = NEXUS ;
fsu.DATAFILE          = Data_file_name   ;
fsu.DATASYN           = DSN   ;
fsu.SYNTHSTYLE        = ST    ;
fsu.SYNTHTREFILE      = STF   ;
fsu.SYNTHTRE          = STR   ;
fsu.TREEPRIOR         = TP    ;
fsu.TOPOLOGYPRIOR     = TOPOLOGYPRIOR; %GKN 18 Mar 2011
fsu.ROOTMAX           = RM    ;
fsu.MCMCINITTREESTYLE = MI    ;
fsu.MCMCINITTREEFILE  = MIF   ;
fsu.MCMCINITTREENUM   = TN    ;
fsu.MCMCINITTREE      = MT    ;
fsu.MCMCINITRHO       = MIR   ;  %RJR 15 Mar 2011
fsu.MCMCINITKAPPA     = MIK   ;  %RJR 15 Mar 2011
fsu.MCMCINITLAMBDA    = 0.1   ;  %ignored anyway
fsu.MCMCCAT           = MCT   ;  %RJR 15 Mar 2011
fsu.MCMCINITMU        = IM    ;  %RJR 15 Mar 2011
fsu.MCMCINITP         = IP     ;  %ignored anyway
fsu.MCMCINITTHETA     = IT    ;
fsu.VARYKAPPA         = MCT   ;  %RJR 15 Mar 2011
fsu.VARYMU            = VARYMU; %RJR 15 Mar 2011
fsu.VARYRHO           = MCT   ;  %RJR 15 Mar 2011
fsu.MCMCMISS          = MISS ; %RJR 15 Mar 2011
fsu.RANDOMKAPPA       = RMIK;
fsu.RANDOMRHO         = RMIR;

fsu.MCMCVARYTOP       = Vary_topology    ;
fsu.VERBOSE           = 1    ; %2 graphs, 1 just textoutput
fsu.OUTFILE           = Output_file_name     ;
fsu.OUTPATH           = Output_path_name    ;
fsu.LOSTONES          = Account_rare_traits    ;
fsu.LOST              = Account_rare_traits    ;
fsu.LOSSRATE          = LR    ;
fsu.ISLRRANDOM        = Random_initial_loss_rate    ;

fsu.PSURVIVE          = PS    ;
fsu.BORROW            = BW    ;
fsu.BORROWFRAC        = BF    ;
fsu.LOCALBORROW       = OFF   ;
fsu.MAXDIST           = 0     ;
fsu.POLYMORPH         = OFF   ;  %TODO XXX FIX CHECKBOXES
fsu.NMEANINGCLASS     = 0     ;
fsu.NU                = 100   ;  %GKN CHANGE Jan 08
fsu.RHO               = 1e-4  ;  %GKN CHANGE Jan 08
fsu.KAPPA             = 0.5   ;  %GKN CHANGE Jan 08
fsu.LAMBDA            = 0.1   ;  %GKN CHANGE Jan 08
fsu.SYNTHCLADES       = 0     ;  %GKN CHANGE Jan 08
fsu.SYNTHCLADESACCURACY = 0   ;  %GKN CHANGE Jan 08
fsu.SYNTHCLADESTIMES  = 0     ;  %GKN CHANGE Jan 08
fsu.NUMSYNTHCLADES    = 0     ;  %GKN CHANGE Jan 08
fsu.MISDAT            = MISS     ;  %RJR 17�Mar 2011
fsu.LOSSRATEBRANCHVAR = 0     ;  %GKN CHANGE Jan 08
fsu.LOSSRATECLASSVAR  = 0     ;  %GKN CHANGE Jan 08

fsu.MASKING           = Omit_taxa    ;
fsu.DATAMASK          = Omit_taxa_list   ;
fsu.ISCOLMASK         = Omit_traits   ;
fsu.COLUMNMASK        = Omit_trait_list    ;
fsu.NUMSEQ            = NS    ;
fsu.VOCABSIZE         = VS    ;
fsu.THETA             = TH    ;
fsu.ISCLADE           = IC    ;
fsu.CLADE             = CL    ;
fsu.CLADEMASK         = CLM   ;
fsu.CLADEAGESMASK     = CLAM  ;
SC=IC;
fsu.STRONGCLADES      = SC    ;
fsu.GUITRUE           = GT    ;
fsu.GUICONTENT        = GC    ;

% Lateral transfer LJK 23/3/20
fsu.BORROWING = BORROWING;
fsu.VARYBETA = VARYBETA;
fsu.MCMCINITBETA = MCMCINITBETA;
fsu.ISBETARANDOM = ISBETARANDOM;

runmcmc(fsu);

end

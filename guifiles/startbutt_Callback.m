function startbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
GlobalValues;

% Clear persistent variables in SDLT code. LUKE 04/09/2016
clear logLkd2_m patternCounts patternMeans

set(h,'UserData',STOPRUN);
ok = 1;
set(handles.statustxt,'String','Initialising');
set([handles.sampledtxt,handles.timegonetxt,handles.timeremtxt],'String','0');

%SET MCMC parameters
% check that MCMC stats are valid
% define RL, SS
vals = str2double(get([handles.runet handles.sampleet],'String'));
RL = vals(1);SS = vals(2);
if ok && (RL < SS)
    disp('Run Length must be greater than Sample Interval')
    ok = 0;
end

% see if we want to seed the random number generator
% define SR, SE
if ok
    if get(handles.seedrandcb,'Value')
        SR = ON; % seeds the rand so runs can be repeated
        SE = str2double(get(handles.seedet,'String'));
    else
        SR = OFF;
        SE = 0;
    end
end

% check that lossrate value is ok
% define LR, IM, RLR
if ok
    RLR = 0; % RLR = is random initial loss rate
    muhandles=[handles.fixmurb handles.specmurb];
    setmu = get(muhandles,'Value');
    setmu = [setmu{:}];
    VARYMU = 1;
    if setmu(1)
        % chosen a fixed mu value
        LR = str2double( get(handles.muvalfixet,'String'));
       VARYMU = 0;
    elseif setmu(2)
        % specified a starting value to vary mu from
        LR = str2double(get(handles.muvalet,'String'));
    else
        %choose random starting mu and vary
        LR = rand;
        RLR = 1;
    end
    IM = DeathRate(LR); % Initial Mu - converted from proportion to modelled mu
end

% Lateral transfer. LUKE 04/09/2016
[BORROWING, VARYBETA, MCMCINITBETA, ISBETARANDOM] = deal(0);
if ok
    % Are we accounting for lateral transfer with SDLT model?
    BORROWING = handles.allowForLateralTransferCB.Value;

    % If so, is beta fixed and what is its initial value?
    if BORROWING
        VARYBETA = handles.varyLateralTransferRateCB.Value;
        ISBETARANDOM = ~handles.initialiseLateralTransferRateAtSpecifiedValueRB.Value;
        if ISBETARANDOM
            MCMCINITBETA = IM * (0.5 + rand);
        else
            MCMCINITBETA = str2double(handles.initialiseLateralTransferRateAtSpecifiedValueEB.String);
        end
    end
end

% check catastrophe parameters
if ok
    RANDOMKAPPA=0; %is initial kappa random
    RANDOMRHO=0;
    if get(handles.includecatscb,'Value') %Include catastrophes
        MCMCCAT=1;

        % kappa
        kappahandles=[handles.fixkapparb, handles.speckapparb];
        setkappa=get(kappahandles,'Value');
        setkappa=[setkappa{:}];
        VARYKAPPA=1;
        if setkappa(1)
            %chosen a fixed kappa value
            MCMCINITKAPPA=str2double(get(handles.kappavalfixet,'String'));
            VARYKAPPA=0;
        elseif setkappa(2)
            %specified a starting value to vary kappa from
            MCMCINITKAPPA=str2double(get(handles.kappavalet,'String'));
        else
            %choose random starting kappa and vary
            MCMCINITKAPPA = 0.25 + 0.75 * rand; % LUKE 02/10/2017
            RANDOMKAPPA=1;
        end

        %rho
        rhohandles=[handles.fixrhorb handles.specrhorb];
        setrho=get(rhohandles,'Value');
        setrho=[setrho{:}];
        VARYRHO=1;
        if setrho(1)
            %chosen a fixed rho value
            MCMCINITRHO=str2double(get(handles.rhovalfixet,'String'));
            VARYRHO=0;
        elseif setrho(2)
            %specified a starting value to vary rho from
            MCMCINITRHO=str2double(get(handles.rhovalet,'String'));
        else
            %choose random starting rho and vary
            [foo k theta]=LogRhoPrior(1); %get the values of k and theta in LogRhoPrior
            MCMCINITRHO=randG(k,1/theta);
            RANDOMRHO=1;
        end
    else % No catastrophes
        MCMCCAT=0;
        MCMCINITKAPPA=0;
        VARYKAPPA=0;
        MCMCINITRHO=0;
        VARYRHO=0;
    end
end

% check missing data parameters
if ok
    MCMCMISS=get(handles.missingdatacb,'Value');
end


% find what kind of initial tree we have
% define MI, IT, MIF, MT, TN
TN = 0;
if ok
    IT =0;MIF = '';MT = [];
    vals = get([handles.randtreerb,handles.spectreerb,handles.truetreerb],'Value');
    val = find([vals{:}]==1);
    if sum([vals{:}]==1)~=1
        disp('More than one intial tree type chosen')
        ok = 0;
    end
    if ok
        switch val
        case 1
            % use Exptree to generate random initial tree
            MI = EXPSTART;
            IT = str2double(get(handles.initthetaet,'String'));
        case 2
            % use a tree stored in a nexus output file to start
            MI = OLDSTART;
            MIF = [handles.oldstart.path handles.oldstart.file];
% TODO values of handles.tree.path etc not loaded in Language:treefilebutt_Callback()
%            MIF = [handles.tree.path handles.tree.file];
            if isempty(MIF)
                fprintf('\nYou need to specify an output tree file from which to start\n')
                ok=0;
            else
                TN  = str2double(get(handles.numtreeet,'String'));
                if ok && TN <= length(handles.tree.output.trees)
                    % make sure that there
                    MT = pop('state');
                    MT.tree = rnextree(handles.tree.output.trees{TN});
                    %TODO this is dangerous - if the columns of
                    %output.stats change meaning then the 4,5 etc below
                    %will be wrong - relevant because btime has been
                    %dropped in writing the output file so these numbers
                    %are not the original output.stats columns
                    MT.mu = handles.tree.output.stats(4,TN);
                    MT.p = handles.tree.output.stats(5,TN);
                    MT.lambda = handles.tree.output.stats(7,TN);
                    if MCMCCAT && ~isempty(handles.tree.output.cattrees)
                        %if you are using catastrophes and an old tree then
                        %there must be a cattree file and you get the old
                        %catastrophes - the following GKN 6/9/09
                        if setkappa(1)
                             MT.kappa=MCMCINITKAPPA;
                        else
                            MT.kappa = handles.tree.output.stats(8,TN);
                        end
                        if setrho(1)
                            MT.rho= MCMCINITRHO;
                        else
                            MT.rho = handles.tree.output.stats(9,TN);
                        end
                        s=MT.tree;
                        sc=rnextree(handles.tree.output.cattrees{TN});
                        %removed legacy code for cattreetolist in svn version < 146
                        [MT.tree MT.cat]=CatTreeToList(sc,s); % RJR 17�Mar 2011
                    else
                        MT.kappa = MCMCINITKAPPA;
                        MT.rho = MCMCINITRHO;
                        MT.cat = [];
                    end
                    % estimate theta at 1/E[edge length]
                    IT = (length(MT.tree)-2)/TreeLength(MT.tree,find([MT.tree.type]==ROOT));

                    % Lateral transfer. LUKE 04/09/2016
                    if BORROWING
                        % We initialise beta = 0 in pop('state').
                        MT.beta = handles.tree.output.stats(12, TN);
                        warning('Catastrophe locations different as not stored in .cat file');
                    end
                else
                    fprintf('\nTree number %1.0f does not exist in the file %s \n',TN,MIF)
                    ok=0;
                end
            end
        case 3
            % use true tree to start from
            MI = TRUSTART;
            IT = handles.data.true.theta;
            MT = pop('state');
            MT.tree = handles.data.true.state.tree;
            MT.mu = handles.data.true.mu;
            MT.p = handles.data.true.p;
            if IT == 0
                % estimate theta at 1/E[edge length]
                IT = (length(MT.tree)-2)/TreeLength(MT.tree,find([MT.tree.type]==ROOT));
            end
            if ~isempty(MT.mu) && LossRate(handles.data.true.mu)>0 && setmu(2)
                % specified a starting value to vary mu from
                IM=MT.mu;
                fprintf('Ignoring the fixed trait death rate set as starting value (%g)\n',LossRate(handles.data.true.mu));
                fprintf('Using the trait death rate %g imported with the true tree to initialise MCMC',LossRate(handles.data.true.mu));
            end

            % Lateral transfer. LUKE 04/09/2016
            MT.beta = handles.data.true.beta;
            if BORROWING && MT.beta > 0
                fprintf('Ignoring the trait transfer rate set as starting value (%g) and using\n', MCMCINITBETA);
                fprintf('the rate %g imported with the true tree to initialise the MCMC instead', MT.beta);
                MCMCINITBETA = MT.beta;
                warning('Catastrophe locations different as not stored in .cat file')
            end
        end
    end
end

% set output parameters
% need to define VB, OP and OF
if ok
    vals = get([handles.drawtreescb handles.plotstatscb handles.quietcb],'Value');
    [dt,ds,qu] = deal(vals{:});
    if qu
        VB = QUIET;
    elseif dt && ds
        VB = GRAPH;
    elseif dt
        VB = JUSTT;
    elseif ds
        VB = JUSTS;
    else
        VB = COUNT;
    end
    % make sure that we have an output file
    if isempty(handles.output.file)
        % use defaults
        OF = get(handles.outfiletxt,'String');
        OP = get(handles.outdirtxt,'String');
    else
        OF = handles.output.file;
        OP = handles.output.path;
    end

end

% find what kind of data we have
% % need to define DS, DSN, DFS, ST, STF, NS, VS, TH, GT, GC
if ok
    if handles.data.truepresent
        DSN = ON;
        GT = handles.data.true;  % GUITRUE is is true state
        GT.NS = length(GT.state.tree)/2;
    else
        DSN = OFF;
        GT = [];
    end
    DS = NEXUS;
    GC = pop('content');  % GUICONTENT is the data from the nexus file
    GC.array = handles.data.array;
    GC.language =handles.data.language;
    GC.cognate = handles.data.cognate;
    [GC.NS,GC.L] = size(GC.array);
    NS = size(handles.data.array,1);
    if isempty(handles.data.file)
        disp('No data file is loaded')
        ok = 0;
% These next 4 lines commented out by RJR 23/05/11. What are they for?
%     elseif isempty(handles.data.array)==0
%         disp(['No data in file ' handles.data.path handles.data.file])
%         DFS = [handles.data.path handles.data.file];
%         %ok = 0;
    else
        DFS = [handles.data.path handles.data.file];
    end
end
if ok
        % called from gui these parameters are arbirtary
    ST = NEWTRE;
    STF = '';
    STR = [];
    VS = 195;
    TH = 1/1000;
    BW = ON;                         %BORROW switch (ON, OFF) simulate borrowing
    BF = 0.15;                       %BORROWFRAC borrowing rate as a fraction of the death rate mu
end

% set model parameters
% need to define TP, RM, LO, LT, VT
if ok
    if get(handles.flatpriorrb,'Value')==1
        TP = FLAT;
    else
        TP = YULE;
    end
    RM = str2double(get(handles.rootet,'String'));
    if get(handles.lostonescb,'Value')
        LO = ON;
    else
        LO = OFF;
    end
    %TODO May reasonably have LT=1 but LO OFF
    LT = LO;
    % are we fixing the tree or not?
    VT = get(handles.varytopcb,'Value');

    if get(handles.uniformpriorrb,'Value')==1
        TOPOLOGYPRIOR=TOPO;
    else
        TOPOLOGYPRIOR=LABHIST;
    end
end


% find out if we are masking any languages
% need to define MK and DM
if ok
    MK = OFF;       %MASK switch {ON OFF}
    DM = [];    %drop these languages from the analysis - its up to you to get the numbers right

    if get(handles.maskcb,'Value') && strcmp(get(handles.maskcb,'Enable'),'on');
        % we are to mask given languages check they are valid
        MK =ON;
        maskstr = get(handles.masket,'String');
        if ~isempty(maskstr)
            mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
            if isempty(mask)
                % str2num could not interpret vector string
                disp([maskstr ' is not a valid vector of taxa to omit'])
                ok = 0;
            else
                % check that all numbers are natural less than number of languages
                if any(floor(mask)~=mask) || mask(1) < 1 || mask(end) > NS
                    fprintf('Taxa to omit must be a vector of integers between 1 and %1.0f \n',NS)
                    ok = 0;
                else
                    DM = mask;
                end
            end
        end
    end
end


%find out if we are masking any cognates
% need to define ICM and CM
if ok
    ICM=OFF; %is clade mask ON or OFF
    CM=[];  % Clade mask
    if get(handles.cogmaskcb,'Value') && strcmp(get(handles.cogmaskcb,'Enable'),'on');
        % we are to mask given languages check they are valid
        ICM =ON;
        maskstr = get(handles.cogmasket,'String');
        if ~isempty(maskstr)
            mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
            if isempty(mask)
                % str2num could not interpret vector string
                disp([maskstr ' is not a valid vector of traits to omit'])
                ok = 0;
            else
                CM = mask;
            end
        end
    end
end


%get clade info
%need to define IC, CL and CLM
if ok
    IC = OFF;
    CL = [];
    CLM = [];
    CLAM=[];
    if (get(handles.cladescb,'Value') == 1) && strcmp(get(handles.cladescb,'Enable'), 'on')
		% impose the clades
        IC = ON;
        CL = handles.data.clade;
        % see whether we need to get rid of any of the clades
        clademaskstr = get(handles.clademasket,'String');
        cladeagesmaskstr=get(handles.cladeagesmasket,'String');
        [ok,IC,CL,CLM,CLAM]=initclades(clademaskstr,IC,CL,cladeagesmaskstr);
    end
end

% make sure that we dont throw out all cognates if we've only 2 languages
if ok && (NS - length(DM)) == 2
    LT = 0;
    LO = OFF;
end
% do some checks that data and initial state match up
%if ok & MK & (MI == TRUSTART)
%    disp('Masking is not implemented when starting from the true tree');
%    ok = 0;
%end
% if ok & IC==ON & (length(CL)>0) & MI == EXPSTART
%     disp(sprintf('\nYou cannot start from a random tree when imposing clades\n'));
%     ok = 0;
% end
if ok && (NS - length(DM)) < 2
     sprintf('\nAt least two taxa must be present for analysis. \n')
     if ~isempty(DM)
		 fprintf('Try shortening the list of omitted taxa.\nCurrently, taxa to be omitted are %s \n',num2str(DM));
     end
     ok = 0;
 end
%if ok & (MI == OLDSTART) & ((MK==OFF) & ((NS - length(DM)) ~= (length(MT.tree)/2)))
%    disp(sprintf('\nData and inital tree don''t match: \n Data has %1.0f languages \n of which %1.0f are to be masked, leaving %1.0f \n but Initial Tree has %1.0f leaves',[NS,length(DM),NS - length(DM),length(MT.tree)/2]));
%    ok = 0;
%end
if ok && (MI == OLDSTART)
    initleafnames = sort({MT.tree(find([MT.tree.type]==LEAF)).Name}'); %#ok<UDIM,FNDSB>
    dataleafnames = sort(handles.data.language(setdiff(1:NS,DM)));
    if (~isequal(initleafnames,dataleafnames) && MK==OFF) || (MK==ON && ~isequal(sort(setdiff(initleafnames,dataleafnames)),sort(handles.data.language(DM))))
        fprintf('\nLeaf names in initial tree don''t match taxon names in data set');
        ok = 0;
    end
end



if ok
    PS=1;
    IP=1;
  %  try
        % run the MCMC unless there is an error

        % configure buttons for run mode
        set(h,'Enable','off');
        set(handles.statustxt,'String','Running');
        set([handles.pausebutt,handles.stopbutt],'Enable','on');

        %write the control variables into structures used by fullsetup
	    fsu=pop('fsu');
        fsu.RUNLENGTH         = RL    ;
        fsu.SUBSAMPLE         = SS    ;
        fsu.SEEDRAND          = SR    ;
        fsu.SEED              = SE    ;
        fsu.DATASOURCE        = DS    ;
        fsu.DATAFILE          = DFS   ;
        fsu.DATASYN           = DSN   ;
        fsu.SYNTHSTYLE        = ST    ;
        fsu.SYNTHTREFILE      = STF   ;
        fsu.SYNTHTRE          = STR   ;
        fsu.TREEPRIOR         = TP    ;
        fsu.ROOTMAX           = RM    ;
        fsu.MCMCINITTREESTYLE = MI    ;
        fsu.MCMCINITTREEFILE  = MIF   ;
        fsu.MCMCINITTREENUM   = TN    ;
        fsu.MCMCINITTREE      = MT    ;
        fsu.MCMCINITMU        = IM    ;
        fsu.MCMCINITP         = IP    ;
        fsu.MCMCINITTHETA     = IT    ;
        fsu.MCMCVARYTOP       = VT    ;
        fsu.VERBOSE           = VB    ;
        fsu.OUTFILE           = OF    ;
        fsu.OUTPATH           = OP    ;
        fsu.LOSTONES          = LO    ;
        fsu.LOST              = LT    ;
        fsu.LOSSRATE          = LR    ;
        fsu.ISLRRANDOM        = RLR    ;
        fsu.PSURVIVE          = PS    ;
        fsu.BORROW            = BW    ;
        fsu.BORROWFRAC        = BF    ;
        fsu.LOCALBORROW       = OFF   ;
        fsu.MAXDIST           = 0     ;
        fsu.POLYMORPH         = OFF   ;  %TODO XXX FIX CHECKBOXES
        fsu.MASKING           = MK    ;
        fsu.DATAMASK          = DM    ;
        fsu.ISCOLMASK         = ICM   ;
        fsu.COLUMNMASK        = CM    ;
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
        fsu.MCMCCAT           = MCMCCAT;
        fsu.RANDOMKAPPA       = RANDOMKAPPA;
        fsu.VARYKAPPA         = VARYKAPPA;
        fsu.MCMCINITKAPPA     = MCMCINITKAPPA;
        fsu.RANDOMRHO         = RANDOMRHO;
        fsu.VARYRHO           = VARYRHO;
        fsu.MCMCINITRHO       = MCMCINITRHO;
        fsu.VARYMU            = VARYMU;
        fsu.MCMCINITLAMBDA    = 1.5e-3; % TODO: remove all instances of lambda everywhere. This value will never be needed (we are integrating lambda out). RJR 19-03-09
        fsu.MCMCMISS          = MCMCMISS;
        MISDAT                = MCMCMISS; % swapped this line and next
        fsu.MISDAT            = MISDAT;   % 19/8/10 GKN
        fsu.TOPOLOGYPRIOR     = TOPOLOGYPRIOR;

        % Lateral transfer. LUKE 04/09/2016
        fsu.BORROWING         = BORROWING;
        fsu.VARYBETA          = VARYBETA;
        fsu.MCMCINITBETA      = MCMCINITBETA;
        fsu.ISBETARANDOM      = ISBETARANDOM;  % LUKE 24/3/20


        runmcmc(fsu,handles,h);
else
    set(handles.statustxt,'String','Idle');
end

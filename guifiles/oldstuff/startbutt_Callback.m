function varargout = startbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
GlobalValues;

set(h,'UserData',STOPRUN); %#ok<NODEF>
ok = 1;
set(handles.statustxt,'String','Initialising');
set([handles.sampledtxt,handles.timegonetxt,handles.timeremtxt],'String','0');

%SET MCMC parameters
% check that MCMC stats are valid
% define RL, SS
vals = str2double(get([handles.runet handles.sampleet],'String'));
RL = vals(1);SS = vals(2);
if ok & (RL < SS)
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
    VARYMU = 1; %#ok<NASGU>
    if setmu(1)
        % chosen a fixed mu value
        LR = str2double( get(handles.muvalfixet,'String'));
       VARYMU = 0; %#ok<NASGU>
    elseif setmu(2)
        % specified a starting value to vary mu from
        LR = str2double(get(handles.muvalet,'String'));
    else
        %chose random starting mu and vary
        LR = rand;
        RLR = 1;
    end
    IM = DeathRate(LR); % Initial Mu - converted from proportion to modelled mu
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
                disp(sprintf('\nYou need to specify an output tree file from which to start'))
                ok=0;
            else
                TN  = str2double(get(handles.numtreeet,'String'));
                if ok & TN <= length(handles.tree.output.trees)
                    % make sure that there 
                    MT = pop('state');
                    MT.tree = rnextree(handles.tree.output.trees{TN});
                    MT.mu = handles.tree.output.stats(4,TN);
                    MT.p = handles.tree.output.stats(5,TN);
                    % estimate theta at 1/E[edge length]
                    IT = (length(MT.tree)-2)/TreeLength(MT.tree,find([MT.tree.type]==ROOT));
                else
                    disp(sprintf('\nTree number %1.0f does not exist in the file %s',TN,MIF))
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
            if ~isempty(MT.mu) & LossRate(handles.data.true.mu)>0 & setmu(2)
                % specified a starting value to vary mu from
                IM=MT.mu;
                disp(sprintf('Ignoring the fixed trait death rate set as starting value (%g)',LossRate(handles.data.true.mu)));
                disp(sprintf('Using the trait death rate %g imported with the true tree to initialise MCMC',LossRate(handles.data.true.mu)));
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
    elseif dt & ds
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
    elseif isempty(handles.data.array)==0
        disp(['No data in file ' handles.data.path handles.data.file])
        DFS = [handles.data.path handles.data.file];
        %ok = 0;
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
end


% find out if we are masking any languages
% need to define MK and DM
if ok
    MK = OFF;       %MASK switch {ON OFF}
    DM = [];    %drop these languages from the analysis - its up to you to get the numbers right

    if get(handles.maskcb,'Value') & strcmp(get(handles.maskcb,'Enable'),'on');
        % we are to mask given languages check they are valid
        MK =ON;
        maskstr = get(handles.masket,'String');
        if ~isempty(maskstr)
            mask = sort(unique(str2num(maskstr)));
            if isempty(mask)
                % str2num could not interpret vector string
                disp([maskstr ' is not a valid vector of taxa to omit'])
                ok = 0;
            else
                % check that all numbers are natural less than number of languages 
                if any(floor(mask)~=mask) | mask(1) < 1 | mask(end) > NS
                    disp(sprintf('Taxa to omit must be a vector of integers between 1 and %1.0f',NS))
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
    if get(handles.cogmaskcb,'Value') & strcmp(get(handles.cogmaskcb,'Enable'),'on');
        % we are to mask given languages check they are valid
        ICM =ON;
        maskstr = get(handles.cogmasket,'String');
        if ~isempty(maskstr)
            mask = sort(unique(str2num(maskstr)));
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
    if (get(handles.cladescb,'Value') == 1) & (get(handles.cladescb,'Enable') == 'on')
        % impose the clades         
        IC = ON;
        CL = handles.data.clade;
        % see whether we need to get rid of any of the clades
        clademaskstr = get(handles.clademasket,'String');
        [ok,IC,CL,CLM]=initclades(clademaskstr,IC,CL);
    end
end

% make sure that we dont throw out all cognates if we've only 2 languages
if ok & (NS - length(DM)) == 2
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
if ok & (NS - length(DM)) < 2
     disp(sprintf('\nAt least two taxa must be left in analysis.  Shorten the list of omitted taxa.\nCurrently, taxa to be omitted are %s \n',num2str(DM)));    
     ok = 0;
 end
%if ok & (MI == OLDSTART) & ((MK==OFF) & ((NS - length(DM)) ~= (length(MT.tree)/2)))
%    disp(sprintf('\nData and inital tree don''t match: \n Data has %1.0f languages \n of which %1.0f are to be masked, leaving %1.0f \n but Initial Tree has %1.0f leaves',[NS,length(DM),NS - length(DM),length(MT.tree)/2]));
%    ok = 0;
%end
if ok & (MI == OLDSTART) 
    initleafnames = sort({MT.tree(find([MT.tree.type]==LEAF)).Name}');
    dataleafnames = sort(handles.data.language(setdiff(1:NS,DM)));
    if (~isequal(initleafnames,dataleafnames) & MK==OFF) | (MK==ON & ~isequal(sort(setdiff(initleafnames,dataleafnames)),sort(handles.data.language(DM))))
        disp(sprintf('\nLeaf names in initial tree don''t match taxon names in data set'));
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
        SC=IC;
        fsu.STRONGCLADES      = SC    ;
        fsu.GUITRUE           = GT    ;
        fsu.GUICONTENT        = GC    ;
	
        
        runmcmc(fsu,handles,h);
else
    set(handles.statustxt,'String','Idle');
end

[SPECIFY DATA SOURCE]
datafilepath = 
datafilename = 

[APPLY MASKS TO DATA AND CLADES]
taxamask =
traitmask =
clademask = 


[SPECIFY INITIAL TREE STATE FOR MCMC]
initialtree = [1 for random, 2 for tree from output file, 3 for true state from data file]
inittreefilename = 
inittreefilepath = 
inittreeindex = 

[SPECIFY TREE PRIOR]
treeprior = [1 for flat, 2 for yule]
maxrootage = [only necessary with flat prior]

estimatelossrate = [0 to leave rate fixed, 1 to estimate]
initlossrate = [between 0 and 1]
userandinitlossrate = 

[PROGRESS REPORTS]
drawcurrentstate = 
drawparametertrace = 
suppressonscreenoutput = 

[SPECIFY MCMC RUN PARAMETERS AND OUTPUT FILE]
outputfilepath = 
outputfilename = 
runlength = 
subsample = 
userandseed = 
seed = 



function varargout = startbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
GlobalValues;

set(h,'UserData',STOPRUN);
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
% define LR, IM
if ok
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
        %chose random starting mu and vary
        LR = rand;
    end
    IM = DeathRate(LR); % Initial Mu - converted from proportion to modelled mu
end

% find what kind of initial tree we have
% define MI, IT, MIF, MT
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
            MIFc = get([handles.treedirtxt handles.treefiletxt],{'String'});           
            MIF = [MIFc{1},MIFc{2}];
% TODO values of handles.tree.path etc not loaded in Language:treefilebutt_Callback()
%            MIF = [handles.tree.path handles.tree.file]; 
            if isempty(MIF)
                disp(sprintf('\nYou need to specify an output tree file from which to start'))
                ok=0;
            else
                treenum  = str2double(get(handles.numtreeet,'String'));
                if ok & treenum <= length(handles.tree.output.trees)
                    % make sure that there 
                    MT = pop('state');
                    MT.tree = rnextree(handles.tree.output.trees{treenum});
                    MT.mu = handles.tree.output.stats(4,treenum);
                    MT.p = handles.tree.output.stats(5,treenum);
                    % estimate theta at 1/E[edge length]
                    IT = (length(MT.tree)-2)/TreeLength(MT.tree,find([MT.tree.type]==ROOT));
                else
                    disp(sprintf('\nTree number %1.0f does not exist in the file %s',treenum,MIF))
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
                disp(sprintf('Ignoring the fixed trait death rate set as starting value',LossRate(handles.data.true.mu)));
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
    elseif length(handles.data.array)==0
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
% need to define TP, RM, LO, LT
if ok
    if get(handles.priorpu,'Value')==1
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
%need to define IC and CL
if ok 
    IC = OFF;
    CL = [];
    if (get(handles.cladescb,'Value') == 1) & (get(handles.cladescb,'Enable') == 'on')
        % impose the clades         
        IC = ON;
        CL = handles.data.clade;
        % see whether we need to get rid of any of the clades
        clademaskstr = get(handles.clademasket,'String');
        if ~isempty(clademaskstr)
            clademask = sort(unique(str2num(clademaskstr)));
            if isempty(clademask)
                % str2num could not interpret vector string
                disp([clademaskstr ' is not a valid vector of clades to omit'])
                ok = 0;
            else
                % check that all numbers are natural and less than number of clades 
                if any(floor(clademask)~=clademask) | clademask(1) < 1 | clademask(end) > length(CL)
                    disp(sprintf('\nClades to omit must be a vector of integers between 1 and %1.0f\n',length(CL)))
                    ok = 0;
                else
                    keepclade = ones(1,length(CL));
                    keepclade(clademask) = 0;
                    disp(sprintf('\nIgnoring the following clades in analysis:'))
                    for i = 1:length(clademask)
                        disp(CL{clademask(i)}.name)
                    end
                    if length(clademask)==length(CL)
                        disp(sprintf('\nAll clades have been dropped, Clading switched off'));
                        CL=[];
                        IC=OFF;
                    else
                        disp('')
                        CL = CL(logical(keepclade));
                    end 
                end
            end
        end
        if ok
            disp(sprintf('\nImposing clades:'))
            for i = 1:length(CL)
                disp(CL{i}.name)
            end
            disp('')
        end
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
        fsu.MCMCINITTREE      = MT    ;
        fsu.MCMCINITMU        = IM    ;
        fsu.MCMCINITP         = IP    ;
        fsu.MCMCINITTHETA     = IT    ;
        fsu.VERBOSE           = VB    ;
        fsu.OUTFILE           = OF    ;
        fsu.OUTPATH           = OP    ;
        fsu.LOSTONES          = LO    ;
        fsu.LOST              = LT    ;
        fsu.LOSSRATE          = LR    ;
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
        SC=IC;
        fsu.STRONGCLADES      = SC    ;
        fsu.GUITRUE           = GT    ;
        fsu.GUICONTENT        = GC    ;
	
        % intialise variables
        [data,model,state,handles.output,mcmc]=fullsetup(fsu);
        save outIC;
        
        set(h,'Interruptible','off');

        if any(handles.output.verbose==[GRAPH JUSTS JUSTT])
            ShowMCMC(model,state,handles.output,data.true);
        end
  
        start=1;
        
        if ( any(handles.output.verbose==[GRAPH JUSTT]) )& ~isempty(data.true.state) 
            draw(data.true.state.tree,handles.output.truefig,LEAF,'true state'); 
        end
         
        disp(sprintf('\n***running MCMC'));
        if ~(handles.output.verbose==QUIET)
            disp(sprintf('(%d,%f)',0,state.loglkd))
        end
        
        finish=floor(mcmc.runlength/mcmc.subsample);
        timestarted = clock;
        
        set(h,'Interruptible','on');
        
        for t=start:finish
 
            %update the Markov chain (mcmc.subsample) steps
            atime=cputime; 
            [state,pa]=Markov(mcmc,model,state); 
            btime=cputime-atime;   
            
            if STOPRUN
                disp('Run halted from stop button')
                STOPRUN = 0;
                set(h,'UserData',STOPRUN);         
                break
            end
            
            if mcmc.gather, save outMC; end   
            
            
            %Write the current state of the Markov chain to output
            NextSamp=handles.output.Nsamp+1;
            stats = [ state.logprior; state.loglkd; state.tree(state.root).time; state.mu; state.p; btime; state.lambda];
            handles.output.stats(:,NextSamp) = stats;
            handles.output.pa(:,NextSamp) = pa;
            handles.output.trees{NextSamp}=wnextree(state.tree,state.root);
            handles.output.Nsamp=NextSamp;
            
            % make output available to GUI
            guidata(gcbf,handles);
            
            %write out progress reports
            set(handles.sampledtxt,'String',sprintf('%4.0f',handles.output.Nsamp));
            timegone = etime(clock,timestarted)/3600;
            set(handles.timegonetxt,'String',sprintf('%4.2f %s',timegone,' hrs'));
            remtime = timegone*(finish+1-(handles.output.Nsamp-1))/(handles.output.Nsamp-1);
            set(handles.timeremtxt,'String',sprintf('%4.2f %s',remtime,' hrs'));
            
            % make sure that the GUI is only interrupted where we want
            set(h,'Interruptible','off');
            %write out progress reports
            if ~(handles.output.verbose==QUIET), disp([sprintf('(%d,%f)',t,state.loglkd),sprintf(' %4.2f',pa')]); end
            if any(handles.output.verbose==[GRAPH,JUSTS,JUSTT]), ShowMCMC(model,state,handles.output,data.true); end
            set(h,'Interruptible','on');
            
            
            %do some routine checking of the state.tree structure
            if TESTSS & check(state,data.true.state)
                check(state,data.true.state)
                disp('Error from check()');
                keyboard;pause;
            end
        end
        
        if mcmc.monitor.on
            mcmc.monitor.data = profile( 'info' );
            profreport( mcmc.monitor.filename, mcmc.monitor.data );
        end
        set(h,'Interruptible','off');
        if ~(handles.output.verbose==QUIET), ShowMCMC(model,state,handles.output,data.true); end
        set(h,'Interruptible','on');
      
%     catch
%         disp(lasterr)
%     end
    save outMC;
    writeoutput(handles.output);
    set(h,'Enable','on');
    set(handles.statustxt,'String','Idle');
    set([handles.pausebutt,handles.stopbutt],'Enable','off');
    disp(['MCMC run finished at ' datestr(clock)]);
else
    set(handles.statustxt,'String','Idle');
end

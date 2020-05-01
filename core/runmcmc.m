function runmcmc(fsu,handles,h)

%global GRAPH JUSTS JUSTT LEAF QUIET ANST STOPRUN TESTSS WRITEXI
GlobalSwitches;
global LOSTONES BORROWING SAVESTATES;

if nargin>0

    fromgui = nargin > 1;

    % intialise variables
    [data,model,state,handles.output,mcmc]=fullsetup(fsu);

    save outIC;

    if fromgui
        set(h,'Interruptible','off');
    end

    if any(handles.output.verbose==[GRAPH JUSTS JUSTT])
        ShowMCMC(model,state,handles.output,data.true);
    end

    start=1;

    if ( any(handles.output.verbose==[GRAPH JUSTT]) ) && ~isempty(data.true.state)
        for i=1:2*state.NS
            if data.true.state.cat(i), data.true.state.tree(i).cat=rand(data.true.state.cat(i),1)*(data.true.state.tree(data.true.state.tree(i).parent).time-data.true.state.tree(i).time)+data.true.state.tree(i).time; else, data.true.state.tree(i).cat=[]; end
        end
        draw(data.true.state.tree,handles.output.truefig,LEAF,'true state');
    end

    % write details of run in par file
    parfilename = [fsu.OUTPATH fsu.OUTFILE '.par'];
    disp(sprintf('Writing parameter file %s\n',parfilename));
    wroteparfile = writerunparfile(parfilename,fsu,fsu.MCMCINITTREENUM);
    if ~wroteparfile
        disp(sprintf('Encountered problem writing parameter file %s.\n Everything else seems ok.',parfilename));
    end

    % write initial state
    handles.output.cattrees{1}=wnexcattree(state.tree,state.root,state.cat);
    writeoutput(handles.output,1);

    if fsu.MCMCMISS && WRITEXI
        writeXIoutput(handles.output,state,1);
    end

    disp(sprintf('\n***running MCMC'));
    if ~(handles.output.verbose==QUIET)
        disp(sprintf('(Sample%5d, loglkd%12.3f)  1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21',0,state.loglkd)) % Luke - added '  21'.
    end

    finish=floor(mcmc.runlength/mcmc.subsample);
    timestarted = clock;

    if fromgui
        set(h,'Interruptible','on');
    end

    % Saving state for later goodness-of-fit testing
    if exist('SAVESTATES', 'var') && ~isempty(SAVESTATES) && SAVESTATES == 1
      [~, ~] = mkdir('saveStates');
      save(sprintf('saveStates%s%s-%05i', filesep, fsu.OUTFILE, 0), 'state');
    end
else

    load outMC
    LOSTONES=model.observe.LostOnes;
    fromgui = 0;
    start=t;

end

if mcmc.gather, lastsave=timestarted; save outMC; end

for t=start:finish

    % Luke 10/02/14
    % Checking to make sure catastrophe locs aren't getting screwed around.
    % We compare the number of catastrophe locations with the number we expect.
    if BORROWING
        for i = 1:(2 * state.NS)
            if state.cat(i) ~= length(state.tree(i).catloc)
                sprintf('Catastrophe mismatch on <pa(%d), %d>', i, i)
                keyboard; pause;
            end
        end
    end

    if mcmc.gather
        if etime(clock,lastsave)>3600, save outMC; lastsave=clock; end
    end

    %update the Markov chain (mcmc.subsample) steps
    atime=cputime;
    ignoreearlywarn= (t<=3); % Ignore warnings which are not alarming when they occur early in the chain. RJR 12/06/11.
    [state,pa]=Markov(mcmc,model,state,ignoreearlywarn);
    btime=cputime-atime;

    if STOPRUN
        disp('Run halted from stop button')
        STOPRUN = 0;
        set(h,'UserData',STOPRUN);
        break
    end

    %Write the current state of the Markov chain to output
    NextSamp=handles.output.Nsamp+1;
    stats = [state.logprior; state.loglkd; state.tree(state.root).time; state.mu; state.p; btime; state.lambda ; state.kappa ; state.rho ; state.ncat ; state.fullloglkd ; state.beta]; % Luke - added ' ; state.beta' to end.
    handles.output.stats(:,NextSamp) = stats;
    handles.output.pa(:,NextSamp) = pa;
    for i=1:2*state.NS  % I think this is superfluous LJK 28/4/20
        if any(state.tree(i).type==[ANST LEAF]), state.tree(i).cat=rand(state.cat(i),1)*(state.tree(state.tree(i).parent).time-state.tree(i).time)+state.tree(i).time; end
    end
    handles.output.trees{NextSamp}=wnextree(state.tree,state.root);
    handles.output.cattrees{NextSamp}=wnexcattree(state.tree,state.root,state.cat);
    handles.output.Nsamp=NextSamp;

    % write to output file
    writeoutput(handles.output,NextSamp);

    % write values of XI
    if fsu.MCMCMISS && WRITEXI
        writeXIoutput(handles.output,state,NextSamp);
    end

    % Saving state for later goodness-of-fit testing
    if exist('SAVESTATES', 'var') && ~isempty(SAVESTATES) && SAVESTATES == 1
      save(sprintf('saveStates%s%s-%05i', filesep, fsu.OUTFILE, t), 'state');
    end

    if fromgui
        % make output available to GUI
        guidata(gcbf,handles);

        %write out progress reports to GUIdes
        set(handles.sampledtxt,'String',sprintf('%4.0f',handles.output.Nsamp));
        timegone = etime(clock,timestarted)/3600;
        set(handles.timegonetxt,'String',sprintf('%4.2f %s',timegone,' hrs'));
        remtime = timegone*(finish+1-(handles.output.Nsamp-1))/(handles.output.Nsamp-1);
        set(handles.timeremtxt,'String',sprintf('%4.2f %s',remtime,' hrs'));

        % make sure that the GUI is only interrupted where we want
        set(h,'Interruptible','off');
    end

    %write out progress reports
    if ~(handles.output.verbose==QUIET), disp([sprintf('(Sample%5d, loglkd%12.3f)',t,state.loglkd),sprintf(' %4.2f',pa')]); end
    if any(handles.output.verbose==[GRAPH,JUSTS,JUSTT]), ShowMCMC(model,state,handles.output,data.true); end

    if fromgui
        set(h,'Interruptible','on');
    end

    % do some routine checking of the state.tree structure
    if TESTSS && check(state,data.true.state)
        format short g
        check(state,data.true.state)
        disp('Error from check()');
        keyboard;pause;
    end                               % LUKE (maybe turn this on again)
end

if mcmc.monitor.on
    mcmc.monitor.data = profile( 'info' );
    profreport( mcmc.monitor.filename, mcmc.monitor.data );
end
if fromgui
    set(h,'Interruptible','off');
end
if ~(handles.output.verbose==QUIET), ShowMCMC(model,state,handles.output,data.true); end
if fromgui
    set(h,'Interruptible','on');
end

save outMC;

%writeoutput(handles.output);

if fromgui
    set(h,'Enable','on');
    set(handles.statustxt,'String','Idle');
    set([handles.pausebutt,handles.stopbutt],'Enable','off');
end

disp(['MCMC run finished at ' datestr(clock)]);

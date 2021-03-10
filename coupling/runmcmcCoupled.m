function runmcmcCoupled(fsu,handles,h)

%global GRAPH JUSTS JUSTT LEAF QUIET ANST STOPRUN TESTSS WRITEXI
GlobalSwitches;
global LOSTONES BORROWING SAVESTATES;

if nargin>0

    fromgui = nargin > 1;

    % intialise variables
    [data,model,state,handles.output,mcmc]=fullsetup(fsu);

    % LJK 23/2/21 not implemented for coupled MCMC
    % save outIC;

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

    [state_x, state_y] = deal(state);
    % Advance each chain by iters draws from prior
    mcmc_prior = mcmc;
    mcmc_prior.subsample = 1e4;
    [state_x, ~] = MarkovPrior(mcmc_prior, model, state_x, 1);
    [state_y, ~] = MarkovPrior(mcmc_prior, model, state_y, 1);

    % Write initial states
    [handles_x, handles_y] = deal(handles);
    handles_x.output.file = [handles_x.output.file, '_x'];
    handles_y.output.file = [handles_y.output.file, '_y'];

    handles_x = write_initial_state(handles_x, state_x, fsu);
    handles_y = write_initial_state(handles_y, state_y, fsu);

    finish=floor(mcmc.runlength/mcmc.subsample);
    timestarted = clock;

    if fromgui
        set(h,'Interruptible','on');
    end

    % Saving state for later goodness-of-fit testing
    if exist('SAVESTATES', 'var') && ~isempty(SAVESTATES) && SAVESTATES == 1
        error('This still needs to be implemented');
        [~, ~] = mkdir('saveStates');
        save(sprintf('saveStates%s%s-%05i', filesep, fsu.OUTFILE, 0), 'state');
    end
else

    error('This still needs to be implemented');
    load outMC
    LOSTONES=model.observe.LostOnes;
    fromgui = 0;
    start=t;

end

% LJK 23/2/21 not implemented for coupled MCMC
% if mcmc.gather, lastsave=timestarted; save outMC; end

% Advance x chain by lag steps
lag_subsample = fsu.COUPLINGLAG / mcmc.subsample;
for t = start:lag_subsample
    borrowing_check(state_x);
    atime = cputime;
    ignoreearlywarn = (t <= 3);
    [state_x, pa_x] = Markov(mcmc, model, state_x, ignoreearlywarn);
    btime = cputime - atime;
    handles_x = write_mcmc_outputs(handles_x, state_x, btime, fsu, t, pa_x, ...
                                   model, data);
end

% Run chain until x reaches finish iterations or coupling, whichever is largest
while t < finish || ~checkCoupling(state_x, state_y)
    t = t + 1;

    borrowing_check(state_x);
    borrowing_check(state_y);

    % LJK 23/2/21 not implemented for coupled MCMC
    % if mcmc.gather
    %     if etime(clock,lastsave)>3600, save outMC; lastsave=clock; end
    % end

    %update the Markov chain (mcmc.subsample) steps
    atime=cputime;
    ignoreearlywarn= (t<=3); % Ignore warnings which are not alarming when they occur early in the chain. RJR 12/06/11.
    [state_x, state_y, pa_x, pa_y] = MarkovCoupled(mcmc, model, state_x, ...
                                                   state_y, ignoreearlywarn);
    btime=cputime-atime;

    if STOPRUN
        disp('Run halted from stop button')
        STOPRUN = 0;
        set(h,'UserData',STOPRUN);
        break
    end

    % Write outputs and update handles
    handles_x = write_mcmc_outputs(handles_x, state_x, btime, fsu, ...
                                   t, pa_x, model, data);
    handles_y = write_mcmc_outputs(handles_y, state_y, btime, fsu, ...
                                   t - lag_subsample, pa_y, model, data);
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

% LJK 23/2/21 not implemented for coupled MCMC
% save outMC;

%writeoutput(handles.output);

if fromgui
    set(h,'Enable','on');
    set(handles.statustxt,'String','Idle');
    set([handles.pausebutt,handles.stopbutt],'Enable','off');
end

disp(['MCMC run finished at ' datestr(clock)]);

end

% Write initial state before loop starts
function [handles] = write_initial_state(handles, state, fsu)
    global WRITEXI QUIET

    % Update handles now that we've sampled from prior
    % TODO: Create subfunction for similar operation in write_mcmc_outputs
    stats = [state.logprior; state.loglkd; state.tree(state.root).time; ...
             state.mu; state.p; 0; state.lambda; state.kappa; state.rho; ...
             state.ncat; state.fullloglkd; state.beta];
    handles.output.stats(:,1) = stats;
    handles.output.trees{1}=wnextree(state.tree,state.root);
    handles.output.cattrees{1}=wnexcattree(state.tree,state.root,state.cat);
    writeoutput(handles.output,1);

    if fsu.MCMCMISS && WRITEXI
        writeXIoutput(handles.output,state,1);
    end

    disp(sprintf('\n***running MCMC'));
    if ~(handles.output.verbose==QUIET)
        disp(sprintf('(Sample%5d, loglkd%12.3f)  1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18   19   20   21',0,state.loglkd)) % Luke - added '  21'.
    end
end

% Borrowing check inside MCMC loop
function [] = borrowing_check(state)
    % Luke 10/02/14
    % Checking to make sure catastrophe locs aren't getting screwed around.
    % We compare the number of catastrophe locations with the number we expect.
    global BORROWING
    if BORROWING
        for i = 1:(2 * state.NS)
            if state.cat(i) ~= length(state.tree(i).catloc)
                sprintf('Catastrophe mismatch on <pa(%d), %d>', i, i)
                keyboard; pause;
            end
        end
    end
end

% Write MCMC outputs
function [handles] = write_mcmc_outputs(handles, state, btime, fsu, t, pa, ...
                                        model, data)
    global ANST LEAF WRITEXI QUIET GRAPH JUSTS JUSTT TESTSS

    %Write the current state of the Markov chain to output
    NextSamp=handles.output.Nsamp+1;
    stats = [state.logprior; state.loglkd; state.tree(state.root).time; state.mu; state.p; btime; state.lambda ; state.kappa ; state.rho ; state.ncat ; state.fullloglkd ; state.beta]; % Luke - added ' ; state.beta' to end.
    handles.output.stats(:,NextSamp) = stats;
    handles.output.pa(:,NextSamp) = pa;
    % % LJK 28/4/20 I think this is superfluous
    % for i=1:2*state.NS
    %     if any(state.tree(i).type==[ANST LEAF]), state.tree(i).cat=rand(state.cat(i),1)*(state.tree(state.tree(i).parent).time-state.tree(i).time)+state.tree(i).time; end
    % end
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
        error('This still needs to be implemented');
      save(sprintf('saveStates%s%s-%05i', filesep, fsu.OUTFILE, t), 'state');
    end

    % LJK 28/4/20 removed for time being
    % if fromgui
    %     % make output available to GUI
    %     guidata(gcbf,handles);
    %
    %     %write out progress reports to GUIdes
    %     set(handles.sampledtxt,'String',sprintf('%4.0f',handles.output.Nsamp));
    %     timegone = etime(clock,timestarted)/3600;
    %     set(handles.timegonetxt,'String',sprintf('%4.2f %s',timegone,' hrs'));
    %     remtime = timegone*(finish+1-(handles.output.Nsamp-1))/(handles.output.Nsamp-1);
    %     set(handles.timeremtxt,'String',sprintf('%4.2f %s',remtime,' hrs'));
    %
    %     % make sure that the GUI is only interrupted where we want
    %     set(h,'Interruptible','off');
    % end

    %write out progress reports
    if ~(handles.output.verbose==QUIET), disp([sprintf('(Sample%5d, loglkd%12.3f)',t,state.loglkd),sprintf(' %4.2f',pa')]); end
    if any(handles.output.verbose==[GRAPH,JUSTS,JUSTT]), ShowMCMC(model,state,handles.output,data.true); end

    % if fromgui
    %     set(h,'Interruptible','on');
    % end

    % do some routine checking of the state.tree structure
    if TESTSS && check(state,data.true.state)
        format short g
        check(state,data.true.state)
        disp('Error from check()');
        keyboard;pause;
    end                               % LUKE (maybe turn this on again)
end

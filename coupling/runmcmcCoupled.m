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
    % Advance each chain by iters mcmc_prior.subsample from prior
    mcmc_prior = mcmc;
    mcmc_prior.subsample = 1e4;
    % Do not change scalar parameters
    mcmc_prior.update.move([8, 15, 17, 21]) = 0;
    mcmc_prior.update.move = mcmc_prior.update.move ./ sum(mcmc_prior.update.move);
    mcmc_prior.update.cmove = cumsum(mcmc_prior.update.move);
    [state_x, ~] = MarkovPrior(mcmc_prior, model, state_x, 1);
    [state_y, ~] = MarkovPrior(mcmc_prior, model, state_y, 1);

    % Write initial states
    [handles_x, handles_y] = deal(handles);
    handles_x.output.file = [handles_x.output.file, '_x'];
    handles_y.output.file = [handles_y.output.file, '_y'];

    handles_x = runmcmcCoupled.writeInitialState(handles_x, state_x, fsu);
    handles_y = runmcmcCoupled.writeInitialState(handles_y, state_y, fsu);

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
    runmcmcCoupled.borrowingCheck(state_x);
    atime = cputime;
    ignoreearlywarn = (t <= 3);
    [state_x, pa_x] = Markov(mcmc, model, state_x, ignoreearlywarn);
    btime = cputime - atime;
    handles_x = runmcmcCoupled.writeMcmcOutputs(handles_x, state_x, btime, ...
                                                fsu, t, pa_x, model, data);
end

% Housekeeping nodes in y to match x
state_y = housekeeping(state_x, state_y);

% Setup proportion accepted output files
runmcmcCoupled.writeProportionAccepted(zeros(mcmc.update.Nmvs, 1), ...
                                       handles_x.output.path, ...
                                       [handles_x.output.file, 'y']);

% Run chain until x reaches finish iterations or coupling, whichever is largest
while t < finish || ~checkCoupling(state_x, state_y)
    t = t + 1;

    runmcmcCoupled.borrowingCheck(state_x);
    runmcmcCoupled.borrowingCheck(state_y);

    % LJK 23/2/21 not implemented for coupled MCMC
    % if mcmc.gather
    %     if etime(clock,lastsave)>3600, save outMC; lastsave=clock; end
    % end

    %update the Markov chain (mcmc.subsample) steps
    atime=cputime;
    ignoreearlywarn= (t<=3); % Ignore warnings which are not alarming when they occur early in the chain. RJR 12/06/11.
    [state_x, state_y, pa_x, pa_y, pa_xy] ...
        = MarkovCoupled(mcmc, model, state_x, state_y, ignoreearlywarn);
    btime=cputime-atime;

    if STOPRUN
        disp('Run halted from stop button')
        STOPRUN = 0;
        set(h,'UserData',STOPRUN);
        break
    end

    % Write outputs and update handles
    handles_x = runmcmcCoupled.writeMcmcOutputs(handles_x, state_x, btime, ...
                                                fsu, t, pa_x, model, data);
    handles_y = runmcmcCoupled.writeMcmcOutputs(handles_y, state_y, btime, ...
                                                fsu, t - lag_subsample, ...
                                                pa_y, model, data);

    % Write proportion of moves accepted
    runmcmcCoupled.writeProportionAccepted(pa_xy, handles_x.output.path, ...
                                           [handles_x.output.file, 'y'], t);
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

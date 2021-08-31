function handles = writeOutputMCMC(...
        handles, state, btime, fsu, t, pa, model, data)

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
    end
end

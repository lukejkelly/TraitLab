% Write initial state before loop starts
function handles = writeInitialState(handles, state, fsu)
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

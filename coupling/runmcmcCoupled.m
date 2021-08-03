function runmcmcCoupled(fsu)

GlobalSwitches;
global EXPSTART

% intialise variables
[data, model, state, handles.output, mcmc] = fullsetup(fsu);

start = 1;

% write details of run in par file
parfilename = [fsu.OUTPATH, fsu.OUTFILE, '.par'];
fprintf('Writing parameter file %s\n\n', parfilename);
wroteparfile = writerunparfile(parfilename, fsu, fsu.MCMCINITTREENUM);
if ~wroteparfile
    fprintf('Encountered problem writing parameter file %s.\n', parfilename);
end

[state_x, state_y] = deal(state);

if mcmc.initial.setup == EXPSTART
    % Advance each chain by sampling from the prior
    [state_x, ~] = MarkovPrior(mcmc, model, state_x, 1);
    [state_y, ~] = MarkovPrior(mcmc, model, state_y, 1);
end

% Write initial states
[handles_x, handles_y] = deal(handles);
handles_x.output.file = [handles_x.output.file, '_x'];
handles_y.output.file = [handles_y.output.file, '_y'];

handles_x = runmcmcCoupled.writeInitialState(handles_x, state_x, fsu);
handles_y = runmcmcCoupled.writeInitialState(handles_y, state_y, fsu);

finish = floor(mcmc.runlength / mcmc.subsample);

% Advance x chain by lag steps
lag_subsample = fsu.COUPLINGLAG / mcmc.subsample;
for t = start:lag_subsample
    atime = cputime;
    ignoreearlywarn = (t <= 3);
    [state_x, pa_x] = Markov(mcmc, model, state_x, ignoreearlywarn);
    btime = cputime - atime;
    handles_x = runmcmcCoupled.writeOutputMCMC(...
        handles_x, state_x, btime, fsu, t, pa_x, model, data);
end

% Housekeeping nodes in y to match x
state_y = housekeeping(state_x, state_y);

% Setup proportion accepted output files
runmcmcCoupled.writeProportionAccepted(zeros(mcmc.update.Nmvs, 1), ...
                                       handles_x.output.path, ...
                                       [handles_x.output.file, 'y']);

% Run coupled kernel until coupling
isCoupled = checkCoupling(state_x, state_y);
if isCoupled
    error('Initial states should not be coupled');
end

while ~isCoupled
    t = t + 1;
    atime = cputime;
    ignoreearlywarn = (t <= 3);
    [state_x, state_y, pa_x, pa_y, pa_xy, isCoupled] = MarkovCoupled(...
        mcmc, model, state_x, state_y, ignoreearlywarn);
    btime = cputime - atime;
    handles_x = runmcmcCoupled.writeOutputMCMC(...
        handles_x, state_x, btime, fsu, t, pa_x, model, data);
    handles_y = runmcmcCoupled.writeOutputMCMC(...
        handles_y, state_y, btime, fsu, t - lag_subsample, pa_y, model, data);
    runmcmcCoupled.writeProportionAccepted(...
        pa_xy, handles_x.output.path, [handles_x.output.file, 'y'], t);
end

runmcmcCoupled.writeCouplingTime(handles_x, t);

% Run marginal to finish if necessary
while t < finish
    t = t + 1;
    atime = cputime;
    ignoreearlywarn = (t <= 3);
    [state_x, pa_x] = Markov(mcmc, model, state_x, ignoreearlywarn);
    btime = cputime - atime;
    handles_x = runmcmcCoupled.writeOutputMCMC(...
        handles_x, state_x, btime, fsu, t, pa_x, model, data);
    runmcmcCoupled.writeProportionAccepted(...
        pa_x, handles_x.output.path, [handles_x.output.file, 'y'], t);
end

disp(['MCMC run finished at ' datestr(clock)]);

end

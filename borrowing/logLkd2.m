function [intLogLik, logLkd, x_R] = logLkd2(state)
% compute integrated and full log-likelihoods of data given tree and parameters

% defining persistent variables
persistent borPars

% converting tree to an array of events and computing the right-left order
% of the leaves
[tEvents, rl] = stype2Events(state.tree);

% creating parameters used in likelihood calculation: if borPars does not
% exist, we create it then update it as the MCMC run progresses
if isempty(borPars)
    borPars = borrowingParameters(tEvents(end).L);
end

% expected pattern frequencies: as x_T scales linearly with lambda, we can
% integrate it out so calculate x_T with lambda = 1
x_T = solveOverTree(tEvents, [1; state.mu; state.beta; state.kappa], borPars);

% to avoid numerical problems
x_T = max(x_T, eps);

% integrated and full log-likelihoods
[intLogLik, logLkd, x_R] = patternCounts(state, rl, x_T, borPars);

end

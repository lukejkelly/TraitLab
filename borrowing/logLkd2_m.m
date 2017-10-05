function [intLogLik, logLkd, x_T] = logLkd2_m( state )
% Function to computed integrated and full log-likelihoods of data given
% state.tree and associated parameters.

% Defining persistent variables.
persistent borPars

% Converting tree to an array of events and computing the right-left order
% of the leaves.
[tEvents, rl] = stype2Events(state.tree);

% Creating parameters used in likelihood calculation. If borPars does not
% exist, we create it then update it as the MCMC run progresses.
if isempty(borPars); borPars = borrowingParameters( tEvents(end).L ); end

% Calculating pattern means. As x_T scales linearly with lambda and we can
% integrate it out, we calculate x_T with lambda = 1.
x_T = solveOverTree(tEvents, [1; state.mu; state.beta; state.kappa], borPars);

% To avoid numerical problems. %% Is this still needed?
x_T = max(x_T, eps);

% Integrated and full log-likelihoods.
[intLogLik, logLkd] = patternCounts(state, rl, x_T, borPars);

end

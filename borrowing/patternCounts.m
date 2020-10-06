function [intLogLik, logLik] = patternCounts(state, rl, x_T, borPars)
% Function which sorts dataset according to topological ordering from
% stype2Events for the current tree then computes pattern counts. Finally,
% the integrated and full log-likelihoods are returned.

% The log-likelihood (ignoring constants) can be written as follows:
%
%   l(n|x) = -\sum_p x_p + \sum_p n_p \log(x_p).
%
% With missing data, this becomes:
%
%   l'(n|x) = -\sum_p x_p + \sum_{p \in P} n_p \log( x_p \xi(p) ) % logLik_o
%             + \sum{q \in Q \ P} n_q \log( x_q \xi(q) )          % logLik_m
%
% Similarly, the integrated log-likelihood (ignoring constants) is:
%
%   h(n|x) = -( \sum_p n_p ) \log( \sum_p y_p ) + \sum_p n_p \log( y_p ),
%
% and with missing data:
%
%   h'(n|x) = -( \sum_{p \in P} n_p ) \log( \sum_{p \in P} y_p )
%                   + \sum_{p \in P} n_p \log( y_p \xi(p) )
%           = -( \sum_{q \in Q\P} n_q ) \log( \sum_{p \in P} y_p )
%                   + \sum_{q \in Q\P} n_q \log( y_q \xi(q) ).

% Note we use the opposite formulation to the TraitLab manual: here, xi is the
% vector of probabilities of observing the true state of a trait at the
% respective leaves

% Setting up data structures, etc.

% First of all we need to account for the fact that the rl (right-to-left) order
% from stype2Events and the order used to count the patterns in the dataset are
% (potentially) not the same.

% MISDAT indicates whether or not we want to model missing data, if any.
global MISDAT LOSTONES

% To avoid counting the patterns each time, we shall define the outputs of
% observedPatternCounts as persistent variables.
persistent obs miss rl_c

if isempty(obs) && isempty(miss) && isempty(rl_c)

    [obs, miss, rl_c] = observedPatternCounts(state);

end

% Number of leaves.
L = state.NS;

% Are there any partially-observed patterns we want to model?
misdat = isstruct(miss) && MISDAT;

% Using the above, we now need to get the pattern counts for the patterns
% defined by the right-left ordering from stype2Events. First of all, we
% shall get the left-right orders, lr_c and lr, then form the permutation
% matrix P such that lr_c * P = lr.

% Left-right orders.
lr_c = fliplr(rl_c);
lr = fliplr(rl);

% Permutation matrix P = [p_{i, j}], where p_{i, j} = 1 if lr_c(i) =
% lr(j), and 0 otherwise.
P = zeros(L);

for j = 1:L

    P(lr_c == lr(j), j) = 1;

end

% Probabilities of observing true values at respective leaves, oriented so
% that row vector corresponds with branch ordering.
xi = [ state.tree(lr).xi ];

% If there is missing data then a registration process removes those patterns
% which do not correspond to at least a single trait being present.
% As a result, we require a correction term in the likelihood calculation.
if misdat

    % We use P_b to find which xi terms are required. For example, the
    % thinning factor for x_{011} is (1 - xi_2)(1 - xi_3).
    C = sum( x_T .* prod( repmat((1 - xi), 2^L - 1, 1) .^ borPars(L).P_b , 2) );
    % C = sum( x_T .* prod( 1 - bsxfun(@times, xi, borPars(L).P_b), 2) );
    % If each entry in xi is 1, that is we always observe the true state, then
    % this reduces to the setting with no missing data setting so C = 0

else

    % If there is no missing data then the registration process has already
    % been accounted for.
    C = 0;

end

% Accounting for correction when traits present in 1 or less leaves are
% discarded.
if LOSTONES

    if misdat

        % For p = 011, say, we subtract xi_2 * (1 - xi_3) + (1 - xi_2) * xi_3
        % off x_011 in the likelihood calculation.
        % TODO: check stability of this calculation when xi is near boundary
        C = C + sum( prod(1 - bsxfun(@times, xi, borPars(L).P_b), 2) ...
                     .* x_T .* sum(bsxfun(@times, xi, borPars(L).P_b) ...
                     ./ (1 - bsxfun(@times, xi, borPars(L).P_b)), 2) );

    else

        C = C + sum(x_T(borPars(L).S == 1));

    end

end

% Contribution to the full log-likelihood from term not involving pattern
% counts: -(\sum_{p \in P} x_p - C) = -(\sum_{r \in R} x_r).
logLik_c = -(sum(x_T) - C) * state.lambda;

% We now compute the contributions to the integrated and full
% log-likelihoods by the fully and partially observed patterns.

% Fully-observed patterns.
if isstruct(obs)

    % Initialising vector of counts.
    n_To = zeros( size(x_T) );

    % Populating the count vector.
    % We use the parallelised C binary-to-decimal function if its compiled
    % version exists. If not, the corresponding Matlab function.
	% n_To( fastBi2De_c_par(obs.pattern * P) ) = obs.count;
	n_To( bi2de(obs.pattern * P, 'left-msb') ) = obs.count;

    % Contribution of observed patterns to integrated and full
    % log-likelihoods.
    if ~misdat % No missing data - ignore xis and correction term.

        % Integrated log-likelihood of fully-observed data (ignoring
        % constants), correction C for LOSTONES registration
        intLogLik_o = sum(n_To .* log(x_T)) - sum(n_To) * log(sum(x_T) - C);

        % Log-likelihood of fully-observed data (ignoring constants).
        logLik_o = sum(n_To .* log(x_T * state.lambda));

    else % Missing data - include xi's in calculation and correct for
        % registration process.

        % In the following, we use the fact that
        % \sum_{q \in Q} x_q \xi(q) = \sum_{p \in P} x_p.

        % Integrated log-likelihood of fully-observed data (ignoring
        % constants).
        intLogLik_o = sum(n_To .* log(x_T * prod(xi))) ...
                      - sum(n_To) * log(sum(x_T) - C);

        % Log-likelihood of fully-observed data (ignoring constants).
        logLik_o = sum(n_To .* log(x_T * prod(xi) * state.lambda));

    end

else

    % No fully-observed patterns so contribution to log-likelihoods is 0.
    intLogLik_o = 0;
    logLik_o = 0;

end

% Partially observed patterns.
if misdat

    % Initialising vector of pattern counts.
    n_Tm = zeros( size(miss) );

    % Initialising vector of pattern means.
    x_Tm = zeros( size(miss) );

    % Populating n_Tm and x_Tm.
    for k = 1:size(miss, 1)

        % Populating vector of pattern counts.
        n_Tm(k) = miss(k).count;

        % Observation probability.
        v = prod( xi.^(miss(k).pattern * P ~= 2) ...
                  .* (1 - xi).^(miss(k).pattern * P == 2) );

        % Missing pattern means.
        % We use parallelised C function for binary-to-decimal operation if it
        % exists. If not, the Matlab function.
		% x_Tm(k) = sum( x_T( fastBi2De_c_par(miss(k).underlying * P) ) ) * v;
        x_Tm(k) = sum( x_T( bi2de(miss(k).underlying * P, 'left-msb') ) ) * v;

    end

    % Integrated log-likelihood of partially-observed patterns (ignoring
    % constants).
    intLogLik_m = sum(n_Tm .* log(x_Tm)) - sum(n_Tm) * log(sum(x_T) - C);

    % Log-likelihood of partially-observed patterns (ignoring constants).
    logLik_m = sum( n_Tm .* log(x_Tm * state.lambda) );

else

    % If no partially-observed patterns then no contribution to
    % log-likelihood.
    intLogLik_m = 0;
    logLik_m = 0;

end

% Integrated log-likelihood.
intLogLik = intLogLik_o + intLogLik_m;

% Log-likelihood.
logLik = logLik_c + logLik_o + logLik_m;

end

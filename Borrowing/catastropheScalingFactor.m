function [logqCat] = catastropheScalingFactor(state, nstate)
% Contribution to log-ratio of proposal densities due to catastrophes when
% the length of the branch they're on changes. The contribution is
%   sum_(i in g') n'_i * log(t'_p(i) - t'_i) - sum_i n_i * log(t_p(i) - t_i).

% Global variables.
global ROOT

% Indices of nodes below root --- those which can carry catastrophes.
sInds = find( [state.tree.type] < ROOT);
nInds = find([nstate.tree.type] < ROOT);

% Scaling term.
logqCat = ...
    sum(nstate.cat(nInds)' .* log([nstate.tree([nstate.tree(nInds).parent]).time] ...
        - [nstate.tree(nInds).time])) ...
  - sum( state.cat(sInds)' .* log([ state.tree([ state.tree(sInds).parent]).time] ...
        - [ state.tree(sInds).time]));

end

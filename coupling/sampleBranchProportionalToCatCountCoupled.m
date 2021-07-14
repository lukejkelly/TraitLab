function [i_x, i_y] = sampleBranchProportionalToCatCountCoupled(state_x, state_y)
    % Sample from a coupling of discrete distributions p and q where for each
    % node i in state_x.tree,
    %     p_i = state_x.cat(i) / state_x.ncat
    % and q likewise
    [rp, ldp] = getDistributionTerms(state_x);
    [rq, ldq] = getDistributionTerms(state_y);
    [i_x, i_y] = maximalCouplingLog(rp, ldp, rq, ldq);
end

function [r, ld] = getDistributionTerms(state)
    r = @() sampleBranchProportionalToCatCount(state);
    ld = @(i) log(state.cat(i)) - log(state.ncat);
end

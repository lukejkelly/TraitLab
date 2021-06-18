function [i_x, i_y] = sampleBranchProportionalToLengthCoupled(state_x, state_y)
    % Sample from a coupling of discrete distributions p and q where for each
    % node i in s = state_x.tree,
    %     p_i = (s(s(i).parent).time - s(i).time) / state_x.length
    % and q likewise
    % Common random numbers didn't work here as even if p_i = q_i, if
    % p_1 + ... + p_{i - 1} > q_1 + ... + q_i then i would never get selected in
    % both states through the core/AddCat inversion sampling scheme
    [rp, ldp] = getDistributionTerms(state_x);
    [rq, ldq] = getDistributionTerms(state_y);
    [i_x, i_y] = maximalCouplingLog(rp, ldp, rq, ldq);
end

function [r, ld] = getDistributionTerms(state)
    r = @() sampleBranchProportionalToLength(state);
    bl = getBranchLengths(state);
    ld = @(i) log(bl(i)) - log(sum(bl));
end

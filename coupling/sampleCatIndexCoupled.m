function [i_x, i_y, j_x, j_y] = sampleCatIndexCoupled(state_x, state_y)
    % Sample from a coupling of discrete distributions p and q where each is
    % uniform over catastrophes on the tree indexed by location
    % k is the location and j its index on branch i's catlocs
    [rp, ldp] = getDistributionTerms(state_x);
    [rq, ldq] = getDistributionTerms(state_y);
    [k_x, k_y] = maximalCouplingLog(rp, ldp, rq, ldq);
    [i_x, j_x] = sampleCatIndexCoupled.getIndex(state_x, k_x);
    [i_y, j_y] = sampleCatIndexCoupled.getIndex(state_y, k_y);
end


function [r, ld] = getDistributionTerms(state)
    r = @() sampleCatIndexCoupled.catlocSample(state);
    ld = @(x) sampleCatIndexCoupled.catlocLogProb(x, state);
end

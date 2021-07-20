function [j_x, j_y] = sampleCatIndexCoupled(state_x, state_y, i_x, i_y)
    % Sample catastrophe indices from coupling of uniform distributions p and q,
    % where p (q) is uniform over catastrophes locations on i_x (i_y)
    [rp, ldp] = getDistributionTerms(state_x, i_x);
    [rq, ldq] = getDistributionTerms(state_y, i_y);
    [k_x, k_y] = maximalCouplingLog(rp, ldp, rq, ldq);
    j_x = sampleCatIndex.getIndex(state_x, i_x, k_x);
    j_y = sampleCatIndex.getIndex(state_y, i_y, k_y);
end

function [r, ld] = getDistributionTerms(state, i)
    r = @() sampleCatIndex.catlocSample(state, i);
    ld = @(x) sampleCatIndex.catlocLogProb(x, state, i);
end

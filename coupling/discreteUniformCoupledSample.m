function [x, y] = discreteUniformCoupledSample(s_x, s_y)
    % Sample (x, y) from coupling of p = Unif(s_x) and q = Unif(s_y)
    [r_p, lp_p] = getDistributionTerms(s_x);
    [r_q, lp_q] = getDistributionTerms(s_y);

    [x, y] = maximalCouplingLog(r_p, lp_p, r_q, lp_q);
end

function [r, lp] = getDistributionTerms(s)
    r = @() discreteUniformSample(s);
    lp = @(z) discreteUniformLogProb(z, s);
end

function [j_x, j_y, k_x, k_y, FAIL_x, FAIL_y] = getWideDestination(i, r, s_x, s_y)
    % Sample from a maximal coupling of discrete uniform distributions on valid
    % destination edge sets
    v_x = BchooseCoupled.getWideDestination.valid(i, r, s_x);
    v_y = BchooseCoupled.getWideDestination.valid(i, r, s_y);

    [r_x, l_x] = getDistributionTerms(v_x);
    [r_y, l_y] = getDistributionTerms(v_y);

    [j_x, j_y] = maximalCouplingLog(r_x, l_x, r_y, l_y);

    [k_x, FAIL_x] = getOutputs(i, s_x, j_x);
    [k_y, FAIL_y] = getOutputs(i, s_y, j_y);
end

function [r, l] = getDistributionTerms(v)
    r = @() BchooseCoupled.getWideDestination.sample(v);
    l = @(x) BchooseCoupled.getWideDestination.logProb(x, v);
end

function [k, FAIL] = getOutputs(i, s, j)
    k = s(j).parent;
    if any([j, k] == s(i).parent)
        FAIL = 1;
    else
        FAIL = 0;
    end
end

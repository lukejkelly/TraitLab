function [i_x, i_y, j_x, j_y] = sampleCatLocIndexCoupled(state_x, state_y)
    % Sample from a coupling of discrete distributions p and q where each is
    % uniform over catastrophe locations on the tree
    % k is the location and j its index on branch i's catlocs
    [rp, ldp] = getDistributionTerms(state_x);
    [rq, ldq] = getDistributionTerms(state_y);
    [k_x, k_y] = maximalCouplingLog(rp, ldp, rq, ldq);
    [i_x, j_x] = getIndex(state_x, k_x);
    [i_y, j_y] = getIndex(state_y, k_y);
end

function [i, j] = getIndex(state, k)
    i = find(cellfun(@(x) ismembertol(k, x), {state.tree.catloc}));
    if length(i) > 1
        error('Only one catastrophe branch should be selected');
    end
    j = find(ismembertol(state.tree(i).catloc, k));
    if length(j) > 1
        error('Only one catastrophe location should be selected');
    end
end

function [r, ld] = getDistributionTerms(state)
    r = @() sampleCatLocIndexCoupled.rCatloc(state);
    ld = @(x) sampleCatLocIndexCoupled.ldCatloc(x, state);
end

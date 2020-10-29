function [nstate_x, nstate_y, logq_x, logq_y, U_x, U_y] ...
     = SchooseCoupledMaximal(state_x, state_y)

    % Housekeeping means root is the same in both x and y, internal nodes
    % labelled with the same set of indices
    i = randsample(state_x.nodes, 1);

    [iT_x, kT_x, jT_x, a_x, b_x] = getNodeTimesAndRanges(i, state_x);
    [iT_y, kT_y, jT_y, a_y, b_y] = getNodeTimesAndRanges(i, state_y);

    if i == state_x.root
        % newage ~ U[(iT + jT) / 2, 2 * iT - jT]
        [newage_x, newage_y] = maximalCouplingUniform(a_x, b_x, a_y, b_y);

        % Sampling densities are (3 / 2) / (iT or newage - jT)
        logq_x = log(iT_x - jT_x) - log(newage_x - jT_x);
        logq_y = log(iT_y - jT_y) - log(newage_y - jT_y);
    else
        % newage ~ U[jT, kT]
        [newage_x, newage_y] = maximalCouplingUniform(jT_x, kT_x, jT_y, kT_y);
        logq_x = 0;
        logq_y = 0;
    end

    [nstate_x, U_x] = Supdate(state_x, i, newage_x);
    [nstate_y, U_y] = Supdate(state_y, i, newage_y);
end

function [iT, kT, jT, a, b] = getNodeTimesAndRanges(i, state)
    % Times of node i, its parent and its oldest child, ranges for sampling
    % new root time
    s = state.tree;
    iT = s(i).time;
    kT = s(s(i).parent).time;
    jT = max(s(s(i).child).time);
    a = (iT + jT) / 2;
    b = 2 * iT - jT;
end

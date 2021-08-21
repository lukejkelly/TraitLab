function [i, newage_x, newage_y, logq_x, logq_y] ...
        = SchooseCoupled(state_x, state_y)

    % Housekeeping means root is the same in both x and y, internal nodes
    % labelled with the same set of indices
    i = state_x.nodes(ceil(rand * (state_x.NS - 1)));

    [iT_x, kT_x, jT_x, a_x, b_x] = SchooseTime.nodeTimesAndRanges(i, state_x);
    [iT_y, kT_y, jT_y, a_y, b_y] = SchooseTime.nodeTimesAndRanges(i, state_y);

    if i == state_x.root
        % newage ~ U[(iT + jT) / 2, 2 * iT - jT]
        [newage_x, newage_y] = maximalCouplingUniform(a_x, b_x, a_y, b_y);

        % Sampling densities are (2 / 3) / (iT or newage - jT)
        logq_x = log(iT_x - jT_x) - log(newage_x - jT_x);
        logq_y = log(iT_y - jT_y) - log(newage_y - jT_y);
    else
        % newage ~ U[jT, kT]
        [newage_x, newage_y] = maximalCouplingUniform(jT_x, kT_x, jT_y, kT_y);
        logq_x = 0;
        logq_y = 0;
    end
end

function [i, newage, logq] = Schoose(state)
    % Sample node i and new age t_i' ~ U(t_j, t_k) if i not root, and if it is
    % then t_i' ~ U[(t_i + t_j) / 2, 2 * t_i - t_j]

    i = state.nodes(ceil(rand * (state.NS - 1)));
    [iT, kT, jT, a, b] = Schoose.nodeTimesAndRanges(i, state);

    if i == state.root
        %this way to update root should be more adaptive to likelihood
        newage = a + rand * (b - a);
        logq = log(iT - jT) - log(newage - jT);
    else
        newage = jT + rand * (kT - jT);
        logq = 0;
    end
end

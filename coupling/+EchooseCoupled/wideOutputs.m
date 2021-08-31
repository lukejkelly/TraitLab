function [i, j, iP, jP, OK] = wideOutputs(state, h)
    % h is linear index into N * N matrix of valid (i, j) pairs
    if isempty(h)
        [i, j, iP, jP] = deal(-1);
        OK = 0;
    else
        N = 2 * state.NS - 1;
        [i, j] = ind2sub([N, N], h);
        iP = state.tree(i).parent;
        jP = state.tree(j).parent;
        OK = 1;
    end
end

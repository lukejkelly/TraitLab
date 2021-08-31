function [catc, ncat, wc, wn] = branchesAndWeights(state, i, newage)
    % Branches to shuffle catastrophes, weights in current and proposed states
    % If i is the root then its weight is automatically 0
    if ismember(i, state.leaves)
        error('Function should be called on leaf nodes');
    end

    % Branches whose catastrophes will be shuffled
    b = [i, state.tree(i).child];

    catc = state.cat(b);
    ncat = sum(state.cat(b));

    % Current branch lengths/weights
    gc = getBranchLengths(state);
    lc = gc(b);
    wc = lc ./ sum(lc);

    % Proposed branch lengths/weights
    nstate = state;
    nstate.tree(i).time = newage;
    gn = getBranchLengths(nstate);
    ln = gn(b);
    wn = ln ./ sum(ln);
end

function v = narrowTargets(state)
    % Get indices of nodes which could move in a narrow edge exchange
    global ROOT
    s = state.tree;
    i = [state.leaves, state.nodes];
    iP = [s(i).parent];
    v = i([s(iP).type] < ROOT);
end

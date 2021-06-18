function x = rCatloc(state)
    % Sample uniformly over catastrophe locations on tree
    i = sampleBranchProportionalToCatCount(state);
    x = state.tree(i).catloc(ceil(state.cat(i) * rand));
end

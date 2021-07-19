function x = catlocSample(state)
    % Sample uniformly over catastrophe locations on tree
    i = sampleBranchByCatCount(state);
    x = state.tree(i).catloc(ceil(state.cat(i) * rand));
end

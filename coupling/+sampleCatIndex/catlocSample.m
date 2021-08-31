function x = catlocSample(state, i)
    % Sample uniformly over catastrophe locations on branch i
    j = sampleCatIndex(state, i);
    x = state.tree(i).catloc(j);
end

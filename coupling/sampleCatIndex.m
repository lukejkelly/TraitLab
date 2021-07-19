function [i, j] = sampleCatIndex(state)
    % Sample catastrophe branch i and location index j uniformly
    k = sampleCatIndexCoupled.catlocSample(state);
    [i, j] = sampleCatIndexCoupled.getIndex(state, k);
end

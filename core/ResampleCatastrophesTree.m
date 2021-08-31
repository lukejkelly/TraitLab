function [nstate, U, OK, logq] = ResampleCatastrophesTree(state)
    % Sample catastrophes from prior, type depends on model and rho
    cat = ResampleCatastrophesTree.samplePriorCounts(state);
    catloc = ResampleCatastrophesTree.samplePriorLocations(cat);
    [nstate, U, OK, logq] = ResampleCatastrophesTree.getOutputs(...
        state, cat, catloc);
end

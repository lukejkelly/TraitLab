function [nstate, U, OK, logq] = ResampleCatastrophes(state)
    % Sample catastrophes from prior, type depends on model and rho
    cat = ResampleCatastrophes.samplePriorCounts(state);
    catloc = ResampleCatastrophes.samplePriorLocations(cat);
    [nstate, U, OK, logq] = ResampleCatastrophes.getOutputs(state, cat, catloc);
end

function [nstate, U, OK, logq] = ResampleCatastrophesBranch(state)
    % Resample catastrophes from the prior on a randomly chosen branch
    global MCMCCAT
    if ~MCMCCAT
        error('Function should only be called when catastrophes switched on');
    end
    i = sampleBranchByLength(state);
    cat = ResampleCatastrophesBranch.samplePriorCount(state, i);
    catloc = ResampleCatastrophesBranch.samplePriorLocations(cat);

    [nstate, U, OK, logq] = ResampleCatastrophesBranch.getOutputs(...
        state, i, cat, catloc);
end

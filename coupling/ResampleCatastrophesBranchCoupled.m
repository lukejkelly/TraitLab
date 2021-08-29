function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = ResampleCatastrophesBranchCoupled(state_x, state_y)
    global MCMCCAT
    if ~MCMCCAT
        error('Function should only be called when catastrophes switched on');
    end

    % Sample target branches
    [i_x, i_y] = sampleBranchByLengthCoupled(state_x, state_y);

    % Sample catastrophe counts
    rp = @() ResampleCatastrophesBranch.samplePriorCount(state_x, i_x);
    lp = @(cat) ResampleCatastrophesBranch.logPriorCount(state_x, i_x, cat);

    rq = @() ResampleCatastrophesBranch.samplePriorCount(state_y, i_y);
    lq = @(cat) ResampleCatastrophesBranch.logPriorCount(state_y, i_y, cat);

    [cat_x, cat_y] = maximalCouplingLog(rp, lp, rq, lq);

    % Sample locations
    catloc = ResampleCatastrophesBranch.samplePriorLocations(max(cat_x, cat_y));

    % Construct outputs
    [nstate_x, U_x, OK_x, logq_x] = ResampleCatastrophesBranch.getOutputs(...
        state_x, i_x, cat_x, catloc);
    [nstate_y, U_y, OK_y, logq_y] = ResampleCatastrophesBranch.getOutputs(...
        state_y, i_y, cat_y, catloc);
end

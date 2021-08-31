function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = ResampleCatastrophesTreeCoupled(state_x, state_y)

    % Couple sampling catastrophes
    rp = @() ResampleCatastrophesTree.samplePriorCounts(state_x);
    lp = @(cat) logProb(cat, state_x);

    rq = @() ResampleCatastrophesTree.samplePriorCounts(state_y);
    lq = @(cat) logProb(cat, state_y);

    [cat_x, cat_y] = maximalCouplingLog(rp, lp, rq, lq);

    % Couple sampling locations
    catloc = ResampleCatastrophesTree.samplePriorLocations(max(cat_x, cat_y));

    [nstate_x, U_x, OK_x, logq_x] ...
        = ResampleCatastrophesTree.getOutputs(state_x, cat_x, catloc);
    [nstate_y, U_y, OK_y, logq_y] ...
        = ResampleCatastrophesTree.getOutputs(state_y, cat_y, catloc);
end

function ld = logProb(cat, state)
    % We are only coupling counts here so use logPriorCounts as the target
    % distribution for coupling, locations are coupled separately, getOutputs
    % uses logPrior (counts and locations) in the ratio of proposals as we do
    % not call catastropheScalingFactor after this move to account for locations
    pstate = state;
    pstate.cat(:) = cat;
    pstate.ncat = sum(cat);
    ld = ResampleCatastrophesTree.logPriorCounts(pstate);
end

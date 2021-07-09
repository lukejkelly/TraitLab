function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = ResampleCatastrophesCoupled(state_x, state_y)

    % Couple sampling catastrophes
    rp = @() ResampleCatastrophes.samplePriorCounts(state_x);
    lp = @(cat) logProb(cat, state_x);

    rq = @() ResampleCatastrophes.samplePriorCounts(state_y);
    lq = @(cat) logProb(cat, state_y);

    [cat_x, cat_y] = maximalCouplingLog(rp, lp, rq, lq);

    % Couple sampling locations
    catloc = ResampleCatastrophes.samplePriorLocations(max(cat_x, cat_y));

    [nstate_x, U_x, OK_x, logq_x] ...
        = ResampleCatastrophes.getOutputs(state_x, cat_x, catloc);
    [nstate_y, U_y, OK_y, logq_y] ...
        = ResampleCatastrophes.getOutputs(state_y, cat_y, catloc);
end

function ld = logProb(cat, state)
    pstate = state;
    pstate.cat(:) = cat;
    pstate.ncat = sum(cat);
    ld = ResampleCatastrophes.logPriorCounts(pstate);
end

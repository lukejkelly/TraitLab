function catastropheScalingFactorCheck(state, nstate, prior)
    % Checking that catastrophe locations are properly accounted for    
    global BORROWING;
    oldBORROWING = BORROWING;
    BORROWING = 0;
    csf = catastropheScalingFactor(state, nstate);
    n1 = - csf + LogPrior(prior, nstate) - LogPrior(prior, state);
    BORROWING = 1;
    n2 = LogPrior(prior, nstate) - LogPrior(prior, state);
    
    if abs(csf) > 1e-10 && abs(n1 - n2) > 1e-10
        keyboard;
    end
    BORROWING = oldBORROWING;
end
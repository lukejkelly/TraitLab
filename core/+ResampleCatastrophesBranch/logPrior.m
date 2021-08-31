function lp = logPrior(state, i)
    % Log-prior on catastrophes on branch (number and locations)
    global BORROWING
    n = state.cat(i);
    lpc = ResampleCatastrophesBranch.logPriorCount(state, i, n);
    if BORROWING
        lpl = ResampleCatastrophesBranch.logPriorLocations(state, i, n);
    else
        lpl = 0;
    end
    lp = lpc + lpl;
end

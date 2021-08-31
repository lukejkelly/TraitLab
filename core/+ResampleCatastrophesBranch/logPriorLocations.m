function lp = logPriorLocations(state, i, n)
    % Log-density (uniform up to permutations) of locations on branch given count
    global BORROWING
    if BORROWING
        d = branchLength(state, i);
        lp = gammaln(n + 1) - n .* log(d);
    else
        lp = 0;
    end
end

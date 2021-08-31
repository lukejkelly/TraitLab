function lpc = logPrior(state)
    % Log-prior on catastrophe number (and locations)
    global VARYRHO BORROWING

    % Contribution from number on tree (if locations are included)
    N = state.ncat;
    D = state.length;
    if VARYRHO
        % Poisson-Gamma(a, b / (b + D))
        [~, a, k] = LogRhoPrior(0);
        b = 1 / k;
        lpn = a * log(b) - gammaln(a) + gammaln(a + N) - (a + N) * log(b + D);
    else
        % Poisson(L * D)
        R = state.rho;
        lpn = N * log(R) - R * D;
    end

    % Integrating out locations (if relevant)
    if BORROWING
        % Not integrating out locations on branches
        lpl = 0;
    else
        % Integrating out locations on branches
        n = state.cat(:)';
        d = getBranchLengths(state);
        i = find(n);
        % Only non-zero catastrophe count branches contribute
        lpl = sum(n(i) .* log(d(i)) - gammaln(n(i) + 1));
    end

    lpc = lpn + lpl;
end

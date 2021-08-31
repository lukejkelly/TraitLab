function lp = logPriorCount(state, i, n)
    % Log-prior on catastrophe number
    global VARYRHO
    d = branchLength(state, i);
    if VARYRHO
        % Poisson-Gamma(a, b / (b + d))
        [~, a, k] = LogRhoPrior(0);
        b = 1 / k;
        p = b / (b + d);
        lp = PoissonGammaLogProb(n, a, p);
    else
        % Poisson
        lp = log(poisspdf(n, d * state.rho));
    end
end

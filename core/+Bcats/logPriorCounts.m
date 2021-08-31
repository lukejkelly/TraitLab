function lp = logPriorCounts(x, d, rho)
    % Log (marginal) prior on catastrophe number on branch of length d
    global VARYRHO

    if VARYRHO
        % x ~ Poisson-Gamma(a, b / (b + d))
        [~, a, k] = LogRhoPrior(0);
        b = 1 / k;
        p = b / (b + d);
        lp = PoissonGammaLogProb(x, a, p);
    else
        % x ~ Poisson(d * rho) internally on log-scale then exponentiated
        r = d * rho;
        lp = log(poisspdf(x, r));
    end
end

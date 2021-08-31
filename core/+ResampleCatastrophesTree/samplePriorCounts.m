function cat = samplePriorCounts(state)
    % Sample catastrophe numbers on branches given lengths
    global VARYRHO

    d = getBranchLengths(state);
    D = state.length;

    if VARYRHO
        % In fact rho is integrated out analytically so Poisson-Gamma counts
        [~, a, k] = LogRhoPrior(0);
        b = 1 / k;

        p = b / (b + D);
        q = d / (b + D);

        cat = MultiPoissonGammaSample(a, p, q);
    else
        % Counts arise from a Poisson(rho) process along branches
        cat = poissrnd(d * state.rho);
    end
end

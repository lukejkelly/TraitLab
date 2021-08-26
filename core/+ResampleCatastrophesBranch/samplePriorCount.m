function cat = samplePriorCount(state, i)
    % Sample catastrophe count on branch i from prior
    global VARYRHO
    d = branchLength(state, i);
    if VARYRHO
        % In fact rho is integrated out analytically so Poisson-Gamma counts
        [~, a, k] = LogRhoPrior(0);
        b = 1 / k;
        p = b / (b + d);
        cat = PoissonGammaSample(a, p);
    else
        % Counts arise from a Poisson(rho) process along branches
        cat = poissrnd(d * state.rho);
    end
end

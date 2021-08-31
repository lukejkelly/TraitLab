function x = samplePriorCounts(d, rho)
    % Sample catastrophe count on a given branch
    global VARYRHO

    if VARYRHO
        % x ~ Poisson-Gamma(a, b / (b + d))
        [~, a, k] = LogRhoPrior(0);
        b = 1 / k;
        p = b / (b + d);
        x = PoissonGammaSample(a, p);
    else
        % x ~ Poisson(d * rho)
        r = d * rho;
        x = poissrnd(r);
    end
end

function p = pg_pdf(x, p)
    [~, a, ~] = LogRhoPrior(0);
    p = exp(PoissonGammaLogProb(x, a, p));
end

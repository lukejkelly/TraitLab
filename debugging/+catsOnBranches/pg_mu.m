function m = pg_mu(p)
    [~, a, ~] = LogRhoPrior(0);
    m = a * (1 - p) / p;
end

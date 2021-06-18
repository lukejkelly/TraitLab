function x = rPoissonGamma(a, p)
    % Sample x from a Poisson-Gamma(a, p) distribution
    %     p(x) = G(a + x) / [G(a) x!] p^a (1 - p)^x
    % so the number of failures until the a'th successs
    % We follow the description in R's stats::rnbinom whereby
    %     x | l ~ Poisson(mean = l) and l ~ Gamma(a, scale = (1 - p) / p)
    l = gamrnd(a, (1 - p) / p);
    x = poissrnd(l);
end

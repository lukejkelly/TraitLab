function ld = ldPoissonGamma(x, a, p)
    % Log-density of x ~ Poisson-Gamma(a, p)
    %     p(x) = G(a + x) / [G(a) x!] p^a (1 - p)^x
    ld = gammaln(a + x) - gammaln(a) - gammaln(x + 1) + a * log(p) ...
         + x * log(1 - p);
end

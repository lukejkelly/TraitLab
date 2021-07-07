function ld = MultiPoissonGammaLogProb(x, a, p, q)
    % Log-mass function of x = (x1, ..., xn) drawn from Multivariate PG
    % generalisation of NegativeMultinomial(a, p, q); that is
    %     x1 + ... + xn = z ~ PG(a, p0)
    %     x | z ~ MN(z, q / (1 - p))
    % x is a vector of the same length as q, a and p are scalars, sum(q) = 1 - p
    ld = gammaln(a + sum(x)) - gammaln(a) - sum(gammaln(x + 1)) ...
         + a * log(p) + sum(x .* log(q));
end

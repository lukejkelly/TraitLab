function x = MultiPoissonGammaSample(a, p, q)
    % Sample x = (x1, ..., xn) from a multivariate PG(a, p) generalisation of
    % the NegativeMultinomial(a, p, q); that is,
    %     x1 + ... + xn = z ~ PG(a, p0)
    %     x | z ~ MN(z, q / (1 - p))
    % x is a vector of the same length as q, a and p are scalars, sum(q) = 1 - p
    z = PoissonGammaSample(a, p);
    x = mnrnd(z, q(:)' / (1 - p));
end

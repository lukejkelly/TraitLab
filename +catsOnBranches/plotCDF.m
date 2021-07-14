function plotCDF(d, n)
    [~, ~, k] = LogRhoPrior(0);
    b = 1 / k;

    [f, x] = ecdf(n);
    x = x(2:end);
    f = f(2:end);
    g = arrayfun(@(t) catsOnBranches.pg_cdf(t, b / (b + d)), 0:max(x));

    yyaxis left
    plot(x, f, 'o', 0:max(x), g, 'x');

    yyaxis right
    plot(x, f - g(:), '*');
end

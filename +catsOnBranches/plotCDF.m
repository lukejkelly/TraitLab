function plotCDF(d, n)
    [~, ~, k] = LogRhoPrior(0);
    b = 1 / k;

    x = 0:max(n);
    f = mean(n(:) <= x);
    g = arrayfun(@(t) mean(catsOnBranches.pg_cdf(t, b ./ (b + d))), x);

    yyaxis left
    plot(x, f, 'o', x, g, 'x');
    xlim([min(x) - 0.1, max(x) + 0.1]);
    ylim([min(min(f), min(g)) * 0.95, max(max(f), max(g)) * 1.05]);

    yyaxis right
    plot(x, f(:) - g(:), '*');
    xlim([min(x) - 0.1, max(x) + 0.1]);
    ylim([min(f(:) - g(:)) * 0.95, max(f(:) - g(:)) * 1.05]);
end

function tests = PoissonGammaTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    % Can only do comparison with NegativeBinomial(a, p) for integer a
    n_i = 3;
    n_j = 1e4;
    for i = 1:n_i
        xObs = nan(1, n_j);
        a = i;
        p = rand;
        for j = 1:n_j
            xObs(j) = rPoissonGamma(a, p);
        end
        subplot(n_i, 1, i);
        [f, x] = ecdf(xObs);
        g = cumsum(arrayfun(@(x) exp(ldPoissonGamma(x, a, p)), 0:max(x)));
        h = nbincdf(x, a, p);
        plot(x, f, 'o', 0:max(x), g, 'x', x, h);
        title(sprintf('PoissonGamma(%d, %g)', a, p));
        ylabel('CDF');
        legend('sample ECDF', 'new', 'exact');
    end
    fprintf('Observed and expected proportions from %g samples\n', n_j);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

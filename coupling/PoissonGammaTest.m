function tests = PoissonGammaTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    % Can only do comparison with NegativeBinomial(a, p) for integer a
    n_i = 3;
    n_j = 1e4;
    for i1 = 1:n_i
        a = i1;
        for i2 = 1:n_i
            p = rand;
            xObs = arrayfun(@(~) PoissonGammaSample(a, p), 1:n_j);
            subplot(n_i, n_i, (i1 - 1) * n_i + i2);
            [f, x] = ecdf(xObs);
            g = cumsum(arrayfun(@(x) exp(PoissonGammaLogProb(x, a, p)), ...
                                0:max(x)));
            h = nbincdf(x, a, p);
            plot(x, f, 'o', 0:max(x), g, 'x', x, h);
            title(sprintf('PoissonGamma(%d, %.3g)', a, p));
            ylabel('CDF');
            legend('sample ECDF', 'new', 'exact');
        end
    end
    fprintf('Observed and expected proportions from %g samples\n', n_j);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function teardownOnce(~)
    close;
end

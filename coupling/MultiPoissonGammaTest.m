function tests = MultiPoissonGammaTest
    tests = functiontests(localfunctions);
end

function simplificationTest(testCase)
    % Compare to PoissonGamma when dimension is 1
    n_i = 3;
    n_j = 1e4;
    for i1 = 1:n_i
        a = -log(rand);
        for i2 = 1:n_i
            p = rand;
            q = 1 - p;
            xObs = arrayfun(@(~) MultiPoissonGammaSample(a, p, q), 1:n_j);
            xExp = arrayfun(@(~) PoissonGammaSample(a, p), 1:n_j);

            [fObs, yObs] = ecdf(xObs);
            [fExp, yExp] = ecdf(xExp);

            subplot(n_i, n_i, (i1 - 1) * n_i + i2);
            plot(yObs, fObs, 'o', yExp, fExp, 'x');
            title(sprintf('PoissonGamma(%.3g, %.3g)', a, p));
            ylabel('CDF');
            legend('Observed', 'Expected');
        end
        z = 0:max([xObs, xExp]);
        zObs = arrayfun(@(x) MultiPoissonGammaLogProb(x, a, p, q), z);
        zExp = arrayfun(@(x) PoissonGammaLogProb(x, a, p), z);
        assertEqual(testCase, zObs, zExp, 'AbsTol', 1e-10);
    end
    fprintf('Observed and expected proportions from %g samples\n', n_j);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function distributionTest(testCase)
    % Compare with NegativeMultinomial(a, p, q) for integer a
    n_i = 3;
    n_j = 1e4;
    for i1 = 1:n_i
        a = i1;
        for i2 = 1:n_i
            p = rand^2;
            m = ceil(rand * 10);
            q = sample_q(p, m);

            [lObs, lExp] = deal(nan(n_j, 1));
            for j = 1:n_j
                x = MultiPoissonGammaSample(a, p, q);
                lObs(j) = MultiPoissonGammaLogProb(x, a, p, q);
                lExp(j) = log(nbinpdf(sum(x), a, p)) ...
                          + log(mnpdf(x, q / (1 - p)));
            end
            assertEqual(testCase, lObs, lExp, 'AbsTol', 1e-10);
        end
    end
end

function q = sample_q(p, n)
    r = -log(rand(1, n));
    q = (1 - p) .* r ./ sum(r);
end

function setupOnce(~)
    clf;
end

function teardownOnce(~)
    close;
end

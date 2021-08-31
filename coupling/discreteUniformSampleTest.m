function tests = discreteUniformSampleTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    N = 12;
    n_i = N / 2;
    n_j = 1e4;
    for i = 1:n_i
        n = 2 * i;
        s = randsample(N, n);
        x = arrayfun(@(~) discreteUniformSample(s), 1:n_j);
        pObs = cumsum(mean(x(:) == 1:N));
        pExp = cumsum(arrayfun(@(x) exp(discreteUniformLogProb(x, s)), 1:N));
        pAct = cumsum(ismember(1:N, s) / n);

        nexttile;
        plot(1:N, pObs, 'o', 1:N, pExp, 'x', 1:N, pAct, '^');
        legend('ECDF', 'CDF1', 'CDF2');
        xlabel('x');
        ylabel('CDF');
        title(sprintf('Support size = %d', n));
    end
    fprintf('Observed and expected proportions from %g samples\n', n_j);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function setupOnce(~)
    clf;
end

function teardownOnce(~)
    close;
end

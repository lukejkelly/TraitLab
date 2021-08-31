function tests = sampleTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    n = 10;
    r = 1:n;

    n_i = 4;
    n_j = 1e4;
    for i = 1:n_i
        u = 1 + floor(rand * n);
        v = randsample(n, u);
        cObs = zeros(size(r));

        for j = 1:n_j
            x = BchooseCoupled.getWideDestination.sample(v);
            cObs(x) = cObs(x) + 1;
        end
        pObs = cObs / n_j;

        lExp = arrayfun(@(x) BchooseCoupled.getWideDestination.logProb(x, v), r);
        pExp = exp(lExp);

        pAct = ismember(r, v) / u;

        nexttile;
        plot(r, pObs, 'o', r, pExp, 'x', r, pAct, '^');
        legend('Observed', 'Expected', 'Truth');
        xlabel('Destination branch');
        ylabel('Proportion sampled');
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

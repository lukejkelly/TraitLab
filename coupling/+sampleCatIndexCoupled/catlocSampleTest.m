function tests = catlocSampleTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    n_i = 4;
    n_j = 1e5;
    for i = 1:n_i
        L = 6 + poissrnd(1);
        state = dummyState(L);

        x = arrayfun(@(~) sampleCatIndexCoupled.catlocSample(state), 1:n_j);

        catlocs = zeros(1, state.ncat);
        k = 1;
        for j = find(state.cat(:)')
            catlocs(k:(k + state.cat(j) - 1)) = state.tree(j).catloc;
            k = k + state.cat(j);
        end
        catlocs = sort(catlocs);
        assertEqual(testCase, length(catlocs), state.ncat);

        fObs = cumsum(mean(x(:) == catlocs));

        pX = @(x) exp(sampleCatIndexCoupled.catlocLogProb(x, state));
        fExp = cumsum(arrayfun(pX, catlocs));

        nexttile;
        plot(catlocs, fObs, 'o', catlocs, fExp, 'x');
        legend('ECDF', 'CDF');
        xlabel('Location');
        ylabel('CDF');
        title('%d leaves', L);
    end
    fprintf('Observed and expected proportions from %g samples\n', n_j);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for c = 1:(1 + poissrnd(3))
        state = AddCat(state);
    end
end


% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 1;
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

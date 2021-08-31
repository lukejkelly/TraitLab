function tests = samplePriorCountTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    % Simulate catastrophes from prior, compare with relevant part of LogPrior
    % We check the support for combinations of 0 and 1 catastrophe
    global VARYRHO
    M = 4;
    N = 1e4;
    for VARYRHO = 0:1
        for m = 1:M
            L = 4 + ceil(rand * 8);
            state = dummyState(L);
            i = sampleBranchByLength(state);

            cats = arrayfun(@(~) ResampleCatastrophesBranch.samplePriorCount(state, i), 1:N);

            nexttile;
            [f, xf] = ecdf(cats);
            xg = 0:max(xf);
            g = cumsum(exp(arrayfun(@(n) ResampleCatastrophesBranch.logPriorCount(state, i, n), xg)));
            plot(xf, f, 'o', xg, g, 'x');
            title(sprintf('VARYRHO = %i : L = %i : i = %i', VARYRHO, L, i));
            xlabel('cats on i');
            ylabel('(E)CDF');
            legend('ecdf', 'cdf');
        end
    end
    fprintf('Observed and expected proportions from %g samples\n', N);
    v1 = input('Are these distributions comparable? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-4));
    state.rho = 2 / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT VARYRHO BORROWING
    BORROWING = 0;
    MCMCCAT = 1;
    testCase.TestData.VARYRHO = VARYRHO;
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    global VARYRHO
    VARYRHO = testCase.TestData.VARYRHO;
    close;
end

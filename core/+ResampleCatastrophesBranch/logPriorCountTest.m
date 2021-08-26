function tests = logPriorCountTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % Simulate catastrophes from prior, compare with relevant part of LogPrior
    global VARYRHO
    M = 10;
    N = 100;
    for VARYRHO = 0:1
        [lpObs, lpExp] = deal(nan(M, N));
        for m = 1:M
            state = dummyState(3 + m);
            for n = 0:N
                i = sampleBranchByLength(state);
                lpObs(m, n + 1) = ResampleCatastrophesBranch.logPriorCount(state, i, n);
                d = branchLength(state, i);
                if VARYRHO
                    [~, a, k] = LogRhoPrior(0);
                    b = 1 / k;
                    lpExp(m, n + 1) = PoissonGammaLogProb(n, a, b / (b + d));
                else
                    lpExp(m, n + 1) = log(poisspdf(n, d * state.rho));
                end
            end
        end
        assertEqual(testCase, lpObs, lpExp, 'AbsTol', 1e-12);
    end
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 5e-4));
    state.rho = 3 / state.length;
    state.cat(:) = ResampleCatastrophesTree.samplePriorCounts(state);
    state.ncat = sum(state.cat);
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global VARYRHO
    testCase.TestData.VARYRHO = VARYRHO;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    global VARYRHO
    VARYRHO = testCase.TestData.VARYRHO;
end

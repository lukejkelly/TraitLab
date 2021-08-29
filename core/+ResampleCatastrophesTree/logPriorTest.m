function tests = logPriorTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % Simulate catastrophes from prior, compare with relevant part of LogPrior
    global YULE LABHIST MCMCCAT BORROWING VARYRHO

    % Minimal prior for LogPrior
    prior.type = YULE;
    prior.topologyprior = LABHIST;

    [~, a, k] = LogRhoPrior(0);
    b = 1 / k;

    n_i = 1e3;
    for L = 5:10
        [lpObs, lpExp] = deal(nan(n_i, 4));
        lpTest =  nan(n_i, 1);
        for i = 1:n_i
            state = dummyState(L);
            for BORROWING = 0:1
                for VARYRHO = 0:1
                    j = 1 + 2 * BORROWING + VARYRHO;
                    lpObs(i, j) = ResampleCatastrophesTree.logPrior(state);
                    % Expected as relevant part of LogPrior
                    MCMCCAT = 0;
                    lpExp0 = LogPrior(prior, state);
                    MCMCCAT = 1;
                    lpExp1 = LogPrior(prior, state);
                    if VARYRHO
                        lpExp1 = lpExp1 + a * log(b) - gammaln(a);
                    end
                    lpExp(i, j) = lpExp1 - lpExp0;
                end
            end
            lpTest(i) = MultiPoissonGammaLogProb(...
                state.cat(:)', ...
                a, ...
                b / (b + state.length), ...
                getBranchLengths(state) / (b + state.length));
        end
        assertEqual(testCase, lpObs, lpExp, 'AbsTol', 1e-12);
        assertEqual(testCase, lpObs(:, 2), lpTest, 'AbsTol', 1e-12);
    end
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 5e-4));
    [~, a, k] = LogRhoPrior(0);
    state.rho = gamrnd(a, k);

    state.cat(:) = ResampleCatastrophesTree.samplePriorCounts(state);
    state.ncat = sum(state.cat);
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT VARYRHO
    MCMCCAT = 1;
    testCase.TestData.VARYRHO = VARYRHO;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    global VARYRHO
    VARYRHO = testCase.TestData.VARYRHO;
end

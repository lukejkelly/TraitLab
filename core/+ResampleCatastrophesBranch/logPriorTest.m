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

    M = 50;
    N = 50;

    for VARYRHO = 0:1
        for BORROWING = 0:1
            [lpObs, lpExp] = deal(nan(M, N));
            for m = 1:M
                L = 3 + m;
                state = dummyState(L);
                for n = 1:N
                    i = sampleBranchByLength(state);
                    cat = ResampleCatastrophesBranch.samplePriorCount(state, i);
                    catloc = ResampleCatastrophesBranch.samplePriorLocations(cat);

                    [nstate, ~, ~, ~] = ResampleCatastrophesBranch.getOutputs(...
                        state, i, cat, catloc);

                    lpObs(n, m) = ResampleCatastrophesBranch.logPrior(nstate, i);

                    % Expected as difference in LogPrior output
                    pstate = nstate;
                    for j = 1:(2 * L)
                        if i ~= j
                            pstate.cat(j) = 0;
                            pstate.tree(j).catloc = [];
                        end
                    end
                    pstate.ncat = cat;
                    MCMCCAT = 0;
                    lpExp0 = LogPrior(prior, pstate);
                    MCMCCAT = 1;
                    lpExp1 = LogPrior(prior, pstate);
                    d = branchLength(pstate, i);
                    if VARYRHO
                        lpExp2 = a * log(b) - gammaln(a) ...
                            + (a + cat) * (log(b + pstate.length) - log(b + d));
                    else
                        lpExp2 = state.rho * (state.length - d);
                    end
                    lpExp(n, m) = lpExp2 + lpExp1 - lpExp0;
                end
            end
            assertEqual(testCase, lpObs, lpExp, 'AbsTol', 1e-12);
        end
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
    global MCMCCAT VARYRHO
    MCMCCAT = 1;
    testCase.TestData.VARYRHO = VARYRHO;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    global VARYRHO
    VARYRHO = testCase.TestData.VARYRHO;
end

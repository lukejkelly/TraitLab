function tests = getOutputsTest
    tests = functiontests(localfunctions);
end

function valuesTest(testCase)
    global BORROWING VARYRHO YULE LABHIST

    % Minimal prior for LogPrior
    prior.type = YULE;
    prior.topologyprior = LABHIST;

    for BORROWING = 0:1
        for VARYRHO = 0:1
            for L = 5:15
                for i = 1:10
                    state = dummyState(L);
                    cat = ResampleCatastrophesTree.samplePriorCounts(state);
                    cat0 = ResampleCatastrophesTree.samplePriorCounts(state);
                    catloc = ResampleCatastrophesTree.samplePriorLocations(...
                                 max(cat, cat0));

                    [nstate, U, OK, logqObs] ...
                        = ResampleCatastrophesTree.getOutputs(state, cat, catloc);

                    assertEqual(testCase, nstate.cat, cat(:));
                    assertEqual(testCase, nstate.ncat, sum(cat));
                    for j = 1:(2 * L)
                        if BORROWING
                            assertSize(testCase, nstate.tree(j).catloc, ...
                                       [cat(j) > 0, cat(j)]);
                        else
                            assertEmpty(testCase, nstate.tree(j).catloc);
                        end
                    end

                    assertEqual(testCase, sort(U), ...
                                sort([state.leaves, state.nodes]));
                    assertEqual(testCase, OK, 1);

                    logqExp = LogPrior(prior, state) ...
                              - LogPrior(prior, nstate);
                    assertEqual(testCase, logqObs, logqExp, 'AbsTol', 1e-10);
                end

            end
        end
    end
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-4));
    state.rho = 3 / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    global VARYRHO MCMCCAT
    unitTests.setupOnce(testCase);
    testCase.TestData.VARYRHO = VARYRHO;
    MCMCCAT = 1;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

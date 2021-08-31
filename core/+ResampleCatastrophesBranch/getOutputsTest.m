function tests = getOutputsTest
    tests = functiontests(localfunctions);
end

function valuesTest(testCase)
    global VARYRHO BORROWING
    M = 10;
    N = 10;
    for BORROWING = 0:1
        for VARYRHO = 0:1
            for m = 1:M
                L = 3 + m;
                state = dummyState(L);
                inds = 1:(2 * L);
                for n = 1:N
                    i = sampleBranchByLength(state);
                    cat = ResampleCatastrophesBranch.samplePriorCount(state, i);
                    catloc = ResampleCatastrophesBranch.samplePriorLocations(2 * cat);
                    [nstate, U, OK, logqObs] = ResampleCatastrophesBranch.getOutputs(...
                        state, i, cat, catloc);

                    assertEqual(testCase, nstate.cat(i), cat);
                    assertEqual(testCase, nstate.cat(inds ~= i), state.cat(inds ~= i));
                    assertEqual(testCase, nstate.ncat, state.ncat + cat - state.cat(i));
                    assertEqual(testCase, nstate.ncat, sum(nstate.cat));
                    for j = inds
                        if BORROWING && nstate.cat(j) > 0
                            if i == j
                                assertEqual(testCase, nstate.tree(j).catloc, ...
                                            sort(catloc(1:cat)));
                            else
                                assertEqual(testCase, nstate.tree(j).catloc, ...
                                            state.tree(j).catloc);
                            end
                        else
                            assertEmpty(testCase, nstate.tree(j).catloc);
                        end
                    end

                    assertEqual(testCase, U, above(i, nstate.tree, nstate.root));
                    assertEqual(testCase, OK, 1);

                    logqExp = ResampleCatastrophesBranch.logPrior(state, i) ...
                              - ResampleCatastrophesBranch.logPrior(nstate, i);
                    assertEqual(testCase, logqObs, logqExp, 'AbsTol', 1e-12);
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

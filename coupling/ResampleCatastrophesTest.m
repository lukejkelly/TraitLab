function tests = ResampleCatastrophesTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % Distribution already checked by subfunctions
    global BORROWING VARYRHO

    for L = 5:10
        for i = 1:10
            for BORROWING = 0:1
                for VARYRHO = 0:1
                    state = dummyState(L);
                    [nstate, U, OK, logq] = ResampleCatastrophes(state);
                    assertEqual(testCase, nstate.ncat, sum(nstate.cat));
                    for j = 1:(2  * nstate.NS)
                        if BORROWING
                            assertEqual(testCase, nstate.cat(j), ...
                                        length(nstate.tree(j).catloc));
                        else
                            assertEmpty(testCase, nstate.tree(j).catloc);
                        end
                    end
                    assertNotEmpty(testCase, U);
                    assertEqual(testCase, OK, 1);
                    assertNotEmpty(testCase, logq);
                end
            end
        end
    end
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 5e-4));
    [~, a, k] = LogRhoPrior(0);
    state.rho = a * k;
    for j = 1:poissrnd(3)
        state = AddCat(state);
    end
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

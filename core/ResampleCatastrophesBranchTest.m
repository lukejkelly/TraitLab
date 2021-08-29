function tests = ResampleCatastrophesBranchTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % Distribution already checked by subfunctions
    global BORROWING VARYRHO
    N = 10;
    for BORROWING = 0:1
        for VARYRHO = 0:1
            for L = 5:10
                state = dummyState(L);
                for n = 1:N
                    [nstate, U, OK, logq] = ResampleCatastrophesBranch(state);
                    ic = find(nstate.cat ~= state.cat);
                    if BORROWING
                        il = find(cellfun(@(x, y) ~(isempty(x) && isempty(y)) && isempty(intersect(x, y)), ...
                                          {nstate.tree.catloc}, {state.tree.catloc}));
                    else
                        il = [];
                    end
                    if isempty(ic)
                        assertEqual(testCase, logq, 0, 'AbsTol', 1e-12);
                        assertEqual(testCase, nstate.cat(U(1)), state.cat(U(1)));
                    elseif BORROWING
                        assertEqual(testCase, ic, il);
                        assertEqual(testCase, ic, U(1));
                    end
                    for j = 1:(2 * nstate.NS)
                        if BORROWING
                            if j == il
                                assertEqual(testCase, nstate.cat(j), ...
                                            length(nstate.tree(j).catloc));
                            else
                                assertEqual(testCase, nstate.tree(j).catloc, ...
                                            state.tree(j).catloc);
                            end
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
    state.rho = 3 / state.length;
    for j = 1:(3 + poissrnd(3))
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

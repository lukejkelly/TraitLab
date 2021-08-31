function tests = BchooseTest
    tests = functiontests(localfunctions);
end

function catastropheTest(testCase)
    % check catastrophe-related outputs only
    global BORROWING MCMCCAT NARROW WIDE
    M = 10;
    N = 10;
    for BORROWING = [0, 1]
        for MCMCCAT = [0, 1]
            for mt = [NARROW, WIDE]
                for m = 1:M
                    L = 5 + ceil(rand * 5);
                    [state, theta, prior] = dummyState(L);
                    for n = 1:N
                        [~, ~, ~, newage, ~, ncat, cat, loc] = Bchoose(...
                            state, mt, theta, prior);

                        if isempty(newage) || ~MCMCCAT
                            assertEmpty(testCase, ncat);
                            assertEmpty(testCase, cat);
                            assertEmpty(testCase, loc);
                        else
                            assertGreaterThanOrEqual(testCase, ncat, 0);
                            assertNotEmpty(testCase, cat);
                            if BORROWING
                                assertNotEmpty(testCase, loc);
                            else
                                assertEmpty(testCase, loc);
                            end
                        end
                    end
                end
            end
        end
    end
end

function [state, theta, prior] = dummyState(L)
    global ROOT
    theta = rand * 1e-2;

    s = ExpTree(L, theta);
    clade = synthclades(s, ceil(rand * L / 2), 2, 1 - rand^3);
    rootmax = (1 + rand) * s([s.type] == ROOT).time;
    prior = unitTests.clade2prior(clade, rootmax);

    state = unitTests.dummyState(s);
    state = UpdateClades(state, [state.leaves, state.nodes], ...
                         size(prior.clade, 2));
    for c = 1:poissrnd(4)
        state = AddCat(state);
    end
    state.rho = state.ncat / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    global VARYRHO
    unitTests.setupOnce(testCase);
    testCase.TestData.VARYRHO = VARYRHO;
    VARYRHO = 1;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

function tests = SchooseCatsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Only checking outputs as distributions checked elsewhere
    global BORROWING
    for BORROWING = 0:1
        for L = 8:12
            state = dummyState(L);
            r = state.root;
            for j = 1:length(state.nodes)
                i = state.nodes(j);
                [~, kT, jT, a, b] = SchooseTime.nodeTimesAndRanges(i, state);
                if i == r
                    newage = a + rand * (b - a);
                else
                    newage = jT + rand * (kT - jT);
                end
                [cat, loc, ~] = SchooseCats(state, i, newage);
                iC = [state.tree(i).child];
                assertEqual(testCase, ...
                            cat.i + cat.j + cat.k, ...
                            sum(state.cat([i, iC])));
                if BORROWING
                    assertEqual(testCase, ...
                                structfun(@length, loc), ...
                                [cat.i; cat.j; cat.k]);
                else
                    assertEmpty(testCase, loc);
                end
            end
        end
    end
end

function setupOnce(testCase)
    global MCMCCAT
    unitTests.setupOnce(testCase);
    MCMCCAT = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for i = 1:(5 + poissrnd(5))
        state = AddCat(state);
    end
end

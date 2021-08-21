function tests = SupdateTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global MCMCCAT BORROWING
    for MCMCCAT = 0:1
        for BORROWING = 0:1
            for L = 8:12
                state = dummyState(L);
                [i, newage, ~, cat, loc] = Schoose(state);
                [nstate, U, TOPOLOGY] = Supdate(state, i, newage, cat, loc);

                assertEqual(testCase, nstate.tree(i).time, newage);
                assertEqual(testCase, ...
                            nstate.length, ...
                            TreeLength(nstate.tree, nstate.root), ...
                            'AbsTol', 1e-10);

                if MCMCCAT
                    bInds = [i, state.tree(i).child];
                    cstate = Supdate.moveCats(state, i, cat, loc);
                    assertEqual(testCase, nstate.cat(bInds), cstate.cat(bInds));
                    if BORROWING
                        assertEqual(testCase, ...
                                    {nstate.tree(bInds).catloc}, ...
                                    {cstate.tree(bInds).catloc});
                    end
                end

                if MCMCCAT
                    UExp = above(state.tree(i).child, state.tree, state.root);
                else
                    UExp = above(i, state.tree, state.root);
                end
                assertEqual(testCase, U, UExp);
                assertEqual(testCase, TOPOLOGY, 0);

            end
        end
    end
end

function setupOnce(testCase)
    unitTests.setupOnce(testCase);
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

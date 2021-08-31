function tests = moveCatsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global BORROWING
    for BORROWING = 0:1
        for L = 8:12
            state = dummyState(L);
            s = state.tree;
            [i, ~, ~, cat, loc] = Schoose(state);

            nstate = Supdate.moveCats(state, i, cat, loc);

            bInds = [i, s(i).child];
            assertEqual(testCase, nstate.cat(bInds)', struct2array(cat));

            cInds = setdiff(1:(2 * state.NS), bInds);
            assertEqual(testCase, nstate.cat(cInds), state.cat(cInds));

            slocs = {state.tree.catloc};
            nlocs = {nstate.tree.catloc};
            if BORROWING
                assertEqual(testCase, ...
                    nlocs(bInds), ...
                    cellfun(@sort, struct2cell(loc), 'UniformOutput', false)');
                assertEqual(testCase, nlocs(cInds), slocs(cInds));
            else
                assertTrue(testCase, all(cellfun(@isempty, slocs)));
                assertTrue(testCase, all(cellfun(@isempty, nlocs)));
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

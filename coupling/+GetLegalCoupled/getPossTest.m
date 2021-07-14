function tests = getPossTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global ROOT
    for L = 5:15
        state = dummyState(L);
        s = state.tree;
        r = state.root;
        for i = find([s.type] < ROOT)
            p_i = s(i).parent;
            c_i = s(i).child;
            s_i = setdiff(s(p_i).child, i);
            n_i = setdiff([p_i, c_i, s_i], r);
            [poss, q] = GetLegalCoupled.getPoss(s, i, r);
            assertEqual(testCase, q, length(n_i));
            assertEqual(testCase, sort(poss), sort(n_i));
            for j = poss
                assertTrue(testCase, ...
                           ismember(i, GetLegalCoupled.getPoss(s, j, r)))
            end
        end
    end
end


function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for c = 1:(1 + poissrnd(3))
        state = AddCat(state);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 0;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

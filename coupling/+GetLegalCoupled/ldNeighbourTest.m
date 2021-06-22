function tests = ldNeighbourTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global ROOT
    for L = 5:15
        state = dummyState(L);
        s = state.tree;
        r = state.root;
        inds = find([s.type] < ROOT);
        for i = inds
            [poss, q] = GetLegalCoupled.getPoss(s, i, r);
            for j = inds
                ld = GetLegalCoupled.ldNeighbour(j, poss, q);
                if ismember(j, poss)
                    assertEqual(testCase, ld, -log(q));
                else
                    assertEqual(testCase, ld, -Inf);
                end
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

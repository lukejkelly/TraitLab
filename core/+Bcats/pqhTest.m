function tests = pqhTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global ROOT
    M = 1e2;
    for m = 1:M
        L = 8 + ceil(rand * 8);
        state = unitTests.dummyState(ExpTree(L, 1e-2));
        v = find([state.tree.type] < ROOT);
        [pObs, pExp, qObs, qExp, hObs, hExp] = deal(nan(size(v)));
        for u = 1:length(v)
            i = v(u);
            [pObs(u), qObs(u), hObs(u)] = Bcats.pqh(state, i);
            pExp(u) = state.tree(i).parent;
            qExp(u) = state.tree(pExp(u)).parent;
            c = state.tree(pExp(u)).child;
            hExp(u) = c(c ~= i);
        end
        assertEqual(testCase, pObs, pExp);
        assertEqual(testCase, qObs, qExp);
        assertEqual(testCase, hObs, hExp);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

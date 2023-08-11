function tests = logpostTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    state = unitTests.dummyState(ExpTree(10, 1e-2));
    for c = 1:(1 + poissrnd(3))
        state = AddCat(state);
    end
    assertEqual(testCase, c, state.ncat);
    rhos = (state.ncat + [-1, 0, 1]) / state.length;
    lpObs = arrayfun(@(x) rho.logpost(x, state), rhos);
    a = 1.5;
    b = 2e-4;
    c = state.ncat;
    d = TreeLength(state.tree, state.root);
    lpExp = log(gampdf(rhos, a + c, 1 / (1 / b + d)));
    assertEqual(testCase, lpObs, lpExp, 'AbsTol', 1e-6);
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

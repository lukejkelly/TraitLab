function tests = logpostTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    state = unitTests.dummyState(ExpTree(10, 1e-2));
    n = 1 + poissrnd(20);
    state.L = n;
    x = 5;
    lambdas = arrayfun(@(~) lambda.sample(state, x), 1:10);
    lpObs = arrayfun(@(y) lambda.logpost(y, state, x), lambdas);
    lpExp = log(gampdf(lambdas, n, 1 / x));
    assertEqual(testCase, lpObs, lpExp, 'AbsTol', 1e-6);
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

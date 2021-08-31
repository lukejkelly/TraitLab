function tests = catlocLogProbTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    M = 10;
    for m = 1:M
        L = 6 + poissrnd(2);
        state = dummyState(L);
        [catlocs, inds] = sampleCatIndex.getCats(state);
        assertEqual(testCase, length(catlocs), state.ncat);
        for n = 1:state.ncat
            i = inds(n, 1);
            lp = @(x) sampleCatIndex.catlocLogProb(x, state, i);
            k = catlocs(n);
            assertEqual(testCase, lp(k), -log(state.cat(i)), 'AbsTol', 1e-12);
            assertEqual(testCase, lp(k - 1e-10), -Inf);
            assertEqual(testCase, lp(k + 1e-10), -Inf);
            assertEqual(testCase, lp(rand), -Inf);
        end
        assertEqual(testCase, lp(-1e-10), -Inf);
        assertEqual(testCase, lp(1 + 1e-10), -Inf);
    end
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for c = 1:(5 + poissrnd(5))
        state = AddCat(state);
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

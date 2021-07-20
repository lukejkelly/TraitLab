function tests = getIndexTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    M = 10;
    for m = 1:M
        L = 6 + poissrnd(2);
        state = dummyState(L);
        [catlocs, indsExp] = sampleCatIndex.getCats(state);

        assertEqual(testCase, length(catlocs), state.ncat);
        indsObs = nan(state.ncat, 1);
        for n = 1:state.ncat
            i = indsExp(n, 1);
            k = catlocs(n);
            indsObs(n) = sampleCatIndex.getIndex(state, i, k);
        end
        assertEqual(testCase, indsObs, indsExp(:, 2));
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

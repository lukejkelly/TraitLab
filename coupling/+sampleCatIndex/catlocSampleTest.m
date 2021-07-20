function tests = catlocSampleTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    M = 4;
    N = 1e3;
    for m = 1:M
        L = 6 + poissrnd(1);
        state = dummyState(L);

        [catlocs, ~] = sampleCatIndex.getCats(state);
        catlocs = sort(catlocs);
        assertEqual(testCase, length(catlocs), state.ncat);

        for i = find(state.cat(:)')
            k = arrayfun(@(~) sampleCatIndex.catlocSample(state, i), 1:N);
            fObs = exp(arrayfun(@(x) sampleCatIndex.catlocLogProb(x, state, i), k));
            fExp = ones(1, N) / length(state.tree(i).catloc);
            assertEqual(testCase, fObs, fExp, 'AbsTol', 1e-10);
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
    BORROWING = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

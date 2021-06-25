function tests = getIndexTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    n_i = 10;
    for i = 1:n_i
        L = 6 + poissrnd(2);
        state = dummyState(L);
        indsExp = nan(state.ncat, 2);
        catlocs = nan(state.ncat, 1);
        k = 1;
        for j = find(state.cat(:)')
            jInds = k:(k + state.cat(j) - 1);
            indsExp(jInds, 1) = j;
            indsExp(jInds, 2) = 1:state.cat(j);
            catlocs(jInds) = state.tree(j).catloc;
            k = k + state.cat(j);
        end
        assertEqual(testCase, length(catlocs), state.ncat);
        indsObs = nan(state.ncat, 2);
        for j = 1:state.ncat
            [indsObs(j, 1), indsObs(j, 2)] ...
                = sampleCatIndexCoupled.getIndex(state, catlocs(j));
        end
        assertEqual(testCase, indsObs(:, 1), indsExp(:, 1));
        assertEqual(testCase, indsObs(:, 2), indsExp(:, 2));
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

function tests = catlocLogProbTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    n_i = 10;
    for i = 1:n_i
        L = 6 + poissrnd(2);
        state = dummyState(L);
        catlocs = zeros(1, state.ncat);
        k = 1;
        for j = find(state.cat(:)')
            catlocs(k:(k + state.cat(j) - 1)) = state.tree(j).catloc;
            k = k + state.cat(j);
        end
        assertEqual(testCase, length(catlocs), state.ncat);
        lp = @sampleCatIndexCoupled.catlocLogProb;
        for j = 1:state.ncat
            c_j = catlocs(j);
            assertEqual(testCase, lp(c_j, state), -log(state.ncat), ...
                        'AbsTol', 1e-12);
            assertEqual(testCase, lp(c_j - 1e-10, state), -Inf);
            assertEqual(testCase, lp(c_j + 1e-10, state), -Inf);
            assertEqual(testCase, lp(rand, state), -Inf);
        end
        assertEqual(testCase, lp(-1e-10, state), -Inf);
        assertEqual(testCase, lp(1 + 1e-10, state), -Inf);
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

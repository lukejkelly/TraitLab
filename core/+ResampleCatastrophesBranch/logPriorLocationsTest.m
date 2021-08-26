function tests = logPriorLocationsTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % Simulate catastrophes from prior, compare with relevant part of LogPrior
    global BORROWING
    M = 10;
    N = 100;
    for BORROWING = 0:1
        [lpObs, lpExp] = deal(nan(M, N));
        for m = 1:M
            state = dummyState(3 + m);
            for n = 0:N
                i = sampleBranchByLength(state);
                lpObs(m, n + 1) = ResampleCatastrophesBranch.logPriorLocations(state, i, n);
                if BORROWING
                    lpExp(m, n + 1) = gammaln(n + 1) - n * log(branchLength(state, i));
                else
                    lpExp(m, n + 1) = 0;
                end
            end
        end
        assertEqual(testCase, lpObs, lpExp, 'AbsTol', 1e-12);
    end
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 5e-4));
    state.rho = 3 / state.length;
    state.cat(:) = ResampleCatastrophesTree.samplePriorCounts(state);
    state.ncat = sum(state.cat);
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function tests = wideOutputsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Compare new and old implementations
    for t = 1:1e2
        L = 8 + ceil(rand * 8);
        state = unitTests.dummyState(ExpTree(L, 1e-3));
        N = 2 * state.NS - 1;
        [iObs1, jObs1, iPObs1, jPObs1, OKObs1] = deal(nan(N));
        [iExp1, jExp1, iPExp1, jPExp1, OKExp1] = deal(nan(N));
        for i = 1:N
            for j = 1:N
                h = i + (j - 1) * N;
                [iObs1(i, j), jObs1(i, j), iPObs1(i, j), jPObs1(i, j), OKObs1(i, j)] ...
                    = EchooseCoupled.wideOutputs(state, h);
                iExp1(i, j) = i;
                jExp1(i, j) = j;
                iPExp1(i, j) = state.tree(i).parent;
                jPExp1(i, j) = state.tree(j).parent;
                OKExp1(i, j) = 1;
            end
            assertEqual(testCase, iObs1, iExp1);
            assertEqual(testCase, jObs1, jExp1);
            assertEqual(testCase, iPObs1, iPExp1);
            assertEqual(testCase, jPObs1, jPExp1);
            assertEqual(testCase, OKObs1, OKExp1);
        end
        [iObs0, jObs0, iPObs0, jPObs0, OKObs0] ...
            = EchooseCoupled.wideOutputs(state, []);
        assertEqual(testCase, iObs0, -1);
        assertEqual(testCase, jObs0, -1);
        assertEqual(testCase, iPObs0, -1);
        assertEqual(testCase, jPObs0, -1);
        assertEqual(testCase, OKObs0, 0);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function tests = narrowTargetsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Compare new and old implementations
    for i = 1:1e2
        L = 8 + ceil(rand * 8);
        state = unitTests.dummyState(ExpTree(L, 1e-4));

        vObs = EchooseCoupled.narrowTargets(state);
        vExp1 = narrowTargetsCheck1(state);
        vExp2 = narrowTargetsCheck2(state);

        assertEqual(testCase, length(vObs), 2 * (L - 2));
        assertEqual(testCase, sort(vObs), sort(vExp1));
        assertEqual(testCase, sort(vObs), sort(vExp2));
    end
end

function v = narrowTargetsCheck1(state)
    iP = state.nodes(state.nodes ~= state.root);
    v = unique([state.tree(iP).child]);
end

function v = narrowTargetsCheck2(state)
    global ADAM
    v = setdiff(1:(2 * state.NS), ...
                [find([state.tree.type] == ADAM), state.root, ...
                    state.tree(state.root).child]);
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

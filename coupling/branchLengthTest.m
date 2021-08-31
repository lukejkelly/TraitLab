function tests = branchLengthTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global ROOT
    for i = 1:10
        L = 5 + ceil(rand * 5);
        s = ExpTree(L, 1e-2);

        state.NS = L;
        state.tree = s;
        state.root = find([s.type] == ROOT);
        state.length = TreeLength(s, state.root);

        blObs = arrayfun(@(i) branchLength(state, i), 1:(2 * state.NS));
        blExp = getBranchLengths(state);

        assertEqual(testCase, blObs, blExp, 'AbsTol', 1e-16);
        assertEqual(testCase, sum(blObs), state.length, 'RelTol', 1e-12);
    end
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

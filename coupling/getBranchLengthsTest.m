function tests = getBranchLengthsTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    global ROOT
    for i = 1:10
        L = 5 + ceil(rand * 5);
        s = ExpTree(L, 1e-2);
        
        state.NS = L;
        state.tree = s;
        state.root = find([s.type] == ROOT);
        state.length = TreeLength(s, state.root);

        blObs = getBranchLengths(state);
        blExp = zeros(size(blObs));
        blExp([s.type] < ROOT) = [s([s([s.type] < ROOT).parent]).time] ...
                - [s([s.type] < ROOT).time];
        
        assertEqual(testCase, blObs, blExp);
        assertEqual(testCase, sum(blObs), state.length, 'RelTol', 1e-12);
    end
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

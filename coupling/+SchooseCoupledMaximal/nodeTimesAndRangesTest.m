function tests = nodeTimesAndRangesTest
    % Unit-testing various maximal coupling functions
    tests = functiontests(localfunctions);
end

function rootTest(testCase)
    global ANST ROOT ADAM
    state.tree = ExpTree(10, 1);
    i = find([state.tree.type] == ROOT);

    [iTObs, kTObs, jTObs, aObs, bObs] ...
        = SchooseCoupledMaximal.nodeTimesAndRanges(i, state);
    assertEqual(testCase, iTObs, state.tree(i).time);
    assertEqual(testCase, kTObs, state.tree([state.tree.type] == ADAM).time);
    assertEqual(testCase, jTObs, ...
                max([state.tree([state.tree.type] == ANST).time]));
    assertEqual(testCase, aObs, mean([iTObs, jTObs]));
    assertEqual(testCase, bObs, 2 * iTObs - jTObs);
    assertEqual(testCase, bObs - aObs, 3 * (iTObs - jTObs) / 2, 'AbsTol', ...
                1e-12);

    [iTExp, kTExp, jTExp, aExp, bExp] = getExpected(i, state.tree);
    assertEqual(testCase, iTObs, iTExp);
    assertEqual(testCase, kTObs, kTExp);
    assertEqual(testCase, jTObs, jTExp);
    assertEqual(testCase, aObs, aExp);
    assertEqual(testCase, bObs, bExp);
end

function anstTest(testCase)
    global ANST
    state.tree = ExpTree(10, 1);
    for i = find([state.tree.type] == ANST)
        [iTObs, kTObs, jTObs, aObs, bObs] ...
            = SchooseCoupledMaximal.nodeTimesAndRanges(i, state);
        assertEqual(testCase, iTObs, state.tree(i).time);
        assertEqual(testCase, kTObs, state.tree(state.tree(i).parent).time);
        assertEqual(testCase, jTObs, ...
                    max([state.tree([state.tree(i).child]).time]));
        assertEqual(testCase, aObs, mean([iTObs, jTObs]));
        assertEqual(testCase, bObs, 2 * iTObs - jTObs);
        assertEqual(testCase, bObs - aObs, 3 * (iTObs - jTObs) / 2, ...
                    'AbsTol', 1e-12);

        [iTExp, kTExp, jTExp, aExp, bExp] = getExpected(i, state.tree);
        assertEqual(testCase, iTObs, iTExp);
        assertEqual(testCase, kTObs, kTExp);
        assertEqual(testCase, jTObs, jTExp);
        assertEqual(testCase, aObs, aExp);
        assertEqual(testCase, bObs, bExp);
    end
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function [iT, kT, jT, a, b] = getExpected(i, s)
    iT = s(i).time;
    k = s(i).parent;
    kT = s(k).time;
    j1 = s(i).child(1);
    j2 = s(i).child(2);
    jT = max(s(j1).time, s(j2).time);
    tau = iT - jT;
    delta_a = -0.5 * tau;
    delta_b = tau;
    a = jT + tau + delta_a;
    b = jT + tau + delta_b;
end

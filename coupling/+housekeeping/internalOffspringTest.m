function tests = internalOffspringTest
    tests = functiontests(localfunctions);
end

function singleTest(testCase)
    global LEAF ANST
    L = 16;
    s = ExpTree(L, 1);
    nodes = find([s.type] == ANST);
    for i = nodes
        rObs = housekeeping.internalOffspring(s, i);
        p = progeny(s, i, LEAF);
        q = intersect(p(1, :), nodes);
        rExp = setdiff(q, i);
        assertEqual(testCase, rObs, rExp);
    end
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

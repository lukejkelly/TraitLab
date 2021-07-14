function tests = getNarrowDestinationTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global LEAF ANST ROOT
    s = ExpTree(10, 1);

    for i = find(ismember([s.type], [LEAF, ANST]))
        [jObs, kObs, FAILObs] ...
            = BchooseCoupled.getNarrowDestination(i, s);
        if s(s(i).parent).type == ROOT
            [jExp, kExp] = deal(-1);
            FAILExp = 1;
        else
            kExp = s(s(i).parent).parent;
            jExp = s(kExp).child(3 - s(s(i).parent).sibling);
            FAILExp = 0;
        end
        assertEqual(testCase, jObs, jExp);
        assertEqual(testCase, kObs, kExp);
        assertEqual(testCase, FAILObs, FAILExp);
    end
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

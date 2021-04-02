function tests = matchSiblingsTest
    tests = functiontests(localfunctions);
end

function orderTest(testCase)
    global ROOT
    for L = randsample(2:15, 1e2, true)
        [s1, s2] = deal(ExpTree(L, 1));
        for i = randsample(find([s1.type] < ROOT), L - 1)
            iP = s2(i).parent;
            iS = s2(iP).child(s2(iP).child ~= i);
            s2(iP).child = fliplr(s2(iP).child);
            [s2([i, iS]).sibling] = deal(s2([iS, i]).sibling);
        end
        assertEqual(testCase, s1, housekeeping.matchSiblings(s1, s2));
        assertEqual(testCase, s2, housekeeping.matchSiblings(s2, s1));
    end
end

function siblingTest(testCase)
    s1 = rnextree('(((1, 2), 3), ((4, 5), 6))');
    s2 = rnextree('((3, (2, 1)), (4, (5, 6)))');

    % adam nodes identical but we need internal node sets to be equal
    for i = 1:6
        j1 = find(cellfun(@(x) strcmp(x, num2str(i)), {s1.Name}));
        if i ~= j1
            s1 = housekeeping.swapNodes(s1, i, j1);
        end
        j2 = find(cellfun(@(x) strcmp(x, num2str(i)), {s2.Name}));
        if i ~= j2
            s2 = housekeeping.swapNodes(s2, i, j2);
        end
    end
    % match nodes with common offspring
    p12 = s1(1).parent;
    s2 = housekeeping.swapNodes(s2, s2(1).parent, p12);

    t2 = housekeeping.matchSiblings(s1, s2);

    for i = 1:3
        assertNotEqual(testCase, s1(i).sibling, s2(i).sibling);
        assertEqual(testCase, s1(i).sibling, t2(i).sibling);
    end
    for i = [4, 6]
        assertEqual(testCase, s1(i).sibling, s2(i).sibling);
        assertEqual(testCase, s1(i).sibling, t2(i).sibling);
    end
    assertNotEqual(testCase, s1(5).sibling, s2(5).sibling);
    assertNotEqual(testCase, s1(5).sibling, t2(5).sibling);

    assertEqual(testCase, s1(p12).child, fliplr(s2(p12).child));
    assertEqual(testCase, s1(p12).child, t2(p12).child);

    p123 = s1(3).parent;
    assertEqual(testCase, s1(p123).child, fliplr(s2(p123).child));
    assertEqual(testCase, s1(p123).child, t2(p123).child);
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function tests = getLeafArrayTest
    tests = functiontests(localfunctions);
end

function permutationTest(testCase)
    global LEAF
    L = 10;
    s = ExpTree(L, 1);

    l1 = find([s.type] == LEAF);
    l2 = l1(randperm(L));
    [~, p2] = sort(l2);

    aObs1 = housekeeping.getLeafArray(s, {s(l1).Name});
    aObs2 = housekeeping.getLeafArray(s, {s(l2).Name});
    assertTrue(testCase, all(aObs1 == aObs2(:, p2), 'all'));
end

function outputTest(testCase)
    % One array started from Adam, which we never normally do
    s = rnextree('(((6, (9, 10)), 5), (1, (2, ((7, 8), (3, 4)))))');
    r = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
    aObs = housekeeping.getLeafArray(s, r);
    aExp = zeros(length(s), length(r));
    aExp(1, :) = 1;
    aExp(2, [5, 6, 9, 10]) = 1;
    aExp(3, [6, 9, 10]) = 1;
    aExp(4, 6) = 1;
    aExp(5, [9, 10]) = 1;
    aExp(6, 9) = 1;
    aExp(7, 10) = 1;
    aExp(8, 5) = 1;
    aExp(9, [1, 2, 3, 4, 7, 8]) = 1;
    aExp(10, 1) = 1;
    aExp(11, [2, 3, 4, 7, 8]) = 1;
    aExp(12, 2) = 1;
    aExp(13, [3, 4, 7, 8]) = 1;
    aExp(14, [7, 8]) = 1;
    aExp(15, 7) = 1;
    aExp(16, 8) = 1;
    aExp(17, [3, 4]) = 1;
    aExp(18, 3) = 1;
    aExp(19, 4) = 1;

    assertTrue(testCase, all(aObs == aExp, 'all'));
end

function recursionTest(testCase)
    for L = 2:15
        recursionTest_(testCase, L);
    end
end

function recursionTest_(testCase, L)
    global LEAF
    s = ExpTree(L, 1);
    r = {s([s.type] == LEAF).Name};
    aObs = housekeeping.getLeafArray(s, r);
    aExp = getLeafArrayWithoutRecursion(s, r);
    assertTrue(testCase, all(aObs == aExp, 'all'));
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function a = getLeafArrayWithoutRecursion(s, r)
    global LEAF ROOT
    a = zeros(length(s), length(r));
    for k = 1:length(s)
        if s(k).type == LEAF
            j = cellfun(@(x) strcmp(x, s(k).Name), r);
            i = k;
            while s(i).type <= ROOT
                a(i, j) = 1;
                i = s(i).parent;
            end
        end
    end
end

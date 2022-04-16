function tests = getLeafArrayTest
    tests = functiontests(localfunctions);
end

function permutationTest(testCase)
    global LEAF ROOT
    L = 100;
    s = ExpTree(L, 1);
    root = find([s.type] == ROOT);

    l1 = find([s.type] == LEAF);
    l2 = l1(randperm(L));
    [~, p2] = sort(l2);

    aObs1 = housekeeping.getLeafArray(s, {s(l1).Name}, root);
    aObs2 = housekeeping.getLeafArray(s, {s(l2).Name}, root);
    assertTrue(testCase, all(aObs1 == aObs2(:, p2), 'all'));

    aOld1 = getLeafArrayOld(s, {s(l1).Name});
    aOld2 = getLeafArrayOld(s, {s(l2).Name});
    assertTrue(testCase, all(aOld1 == aOld2(:, p2), 'all'));
    assertTrue(testCase, all(aObs1 == aOld1, 'all'));
    assertTrue(testCase, all(aObs2 == aOld2, 'all'));
end

function outputTest(testCase)
    global ROOT
    s = rnextree('(((6, (9, 10)), 5), (1, (2, ((7, 8), (3, 4)))))');
    r = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'};
    root = find([s.type] == ROOT);
    aObs = housekeeping.getLeafArray(s, r, root);
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

    aOld = getLeafArrayOld(s, r);
    assertTrue(testCase, all(aObs == aOld, 'all'));
end

function recursionTest(testCase)
    for L = 2:50
        recursionTest_(testCase, L);
    end
end

function recursionTest_(testCase, L)
    global LEAF ROOT
    s = ExpTree(L, 1);
    r = {s([s.type] == LEAF).Name};
    root = find([s.type] == ROOT);
    aObs = housekeeping.getLeafArray(s, r, root);
    aExp1 = getLeafArrayWithoutRecursion(s, r);
    assertTrue(testCase, all(aObs == aExp1, 'all'));
    aExp2 = getLeafArrayOld(s, r);
    assertTrue(testCase, all(aObs == aExp2, 'all'));
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

function a = getLeafArrayOld(s, r, a, l)
    % Starting at node l (root or lower) and working towards the leaves of s,
    % a indicates the leaves (identified by name) beneath each node according to
    % the reference list r; that is, a is a |s|x|a| array with a_ij = 1 if leaf
    % r_j is beneath node s_i

    global LEAF ROOT;

    if nargin == 2
        a = zeros(length(s), length(r));
        l = find([s.type] == ROOT);
    end

    if s(l).type == LEAF
        a(l, ismember(r, s(l).Name)) = 1;
    else
        a1 = getLeafArrayOld(s, r, a, s(l).child(1));
        a2 = getLeafArrayOld(s, r, a, s(l).child(2));
        a(~all(a1 == 0, 2), :) = a1(~all(a1 == 0, 2), :);
        a(~all(a2 == 0, 2), :) = a2(~all(a2 == 0, 2), :);
        a(l, :) = sum(a(s(l).child, :));
    end
end

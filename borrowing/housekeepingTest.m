function tests = housekeepingTest
    % Unit-testing likelihood calculations in patternCounts
    tests = functiontests(localfunctions);
end

%% Test Functions
function s1s1Test(testCase)
    % Identical inputs yield identical outputs
    s1 = testCase.TestData.s1;
    t1Obs = housekeeping(s1, s1);
    checkHousekeeping(testCase, t1Obs, s1);
end

function s1s2Test(testCase)
    % Switch root and internal node, overlapping subtree sets
    s1 = testCase.TestData.s1;
    s2 = testCase.TestData.s2;
    t2Obs = housekeeping(s1, s2);
    t2Exp = testCase.TestData.t2s1;
    checkHousekeeping(testCase, t2Obs, t2Exp);
end

function s1s3Test(testCase)
    % Switch root and internal node, non-overlapping subtree sets
    s1 = testCase.TestData.s1;
    s3 = testCase.TestData.s3;
    t3Obs = housekeeping(s1, s3);
    t3Exp = testCase.TestData.t3s1;
    checkHousekeeping(testCase, t3Obs, t3Exp);
end

function s1s4Test(testCase)
    % Switch Adam, root internal node and leaves, non-overlapping subtree sets
    s1 = testCase.TestData.s1;
    s4 = testCase.TestData.s4;
    t4Obs = housekeeping(s1, s4);
    t4Exp = testCase.TestData.t4s1;
    checkHousekeeping(testCase, t4Obs, t4Exp);
end

function s1s5Test(testCase)
    % Permutation of indices, leaves at different times
    s1 = testCase.TestData.s1;
    s5 = testCase.TestData.s5;
    t5Obs = housekeeping(s1, s5);
    t5Exp = testCase.TestData.t5s1;
    checkHousekeeping(testCase, t5Obs, t5Exp);
end

function s6s7Test(testCase)
    % Permutation of indices, leaves at different times
    s6 = testCase.TestData.s6;
    s7 = testCase.TestData.s7;
    t7Obs = housekeeping(s6, s7);
    t7Exp = testCase.TestData.t7s6;
    checkHousekeeping(testCase, t7Obs, t7Exp);
end

function s6s8Test(testCase)
    % Only one node swap, leaves at different times
    s6 = testCase.TestData.s6;
    s8 = testCase.TestData.s8;
    t8Obs = housekeeping(s6, s8);
    t8Exp = testCase.TestData.t8s6;
    checkHousekeeping(testCase, t8Obs, t8Exp);
end

% Helper functions
function checkHousekeeping(testCase, tObs, tExp)
    fn = {'Name', 'parent', 'child', 'type', 'time', 'sibling'};
    for j = 1:length(fn)
        assertTrue(testCase, ...
                   all(cellfun(@(a, b) isequaln(a, b), ...
                               {tObs.(fn{j})}, ...
                               {tExp.(fn{j})})));
    end
    assertTrue(testCase, equaltrees(tObs, tExp) == 1);
end

% Optional file fixtures
function setupOnce(testCase)
    % Setting up
    GlobalSwitches;
    GlobalValues;
    addpath('debugging', '.');

    % Root and left child are sibling 1, right child is sibling 2
    % Leaves are at time 0, internal node at 1, root at 2 and Adam at 3

    % s1
    % 6
    % |
    % 5 — 3
    % |
    % 4 — 2
    % |
    % 1
    testCase.TestData.s1 = struct(...
        'Name', {'1', '2', '3', '', '', 'Adam'}, ...
        'parent', {4, 4, 5, 5, 6, [],}, ...
        'child', {[], [], [], [1, 2], [4, 3], 5}, ...
        'type', {0, 0, 0, 1, 2, 3}, ...
        'time', {0, 0, 0, 1, 2, 3}, ...
        'sibling', {1, 2, 2, 1, 1, []}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % s2
    % 6
    % |
    % 4 — 3
    % |
    % 5 — 2
    % |
    % 1
    testCase.TestData.s2 = struct(...
        'Name', {'1', '2', '3', '', '', 'Adam'}, ...
        'parent', {5, 5, 4, 6, 4, [],}, ...
        'child', {[], [], [], [5, 3], [1, 2], 4}, ...
        'type', {0, 0, 0, 2, 1, 3}, ...
        'time', {0, 0, 0, 2, 1, 3}, ...
        'sibling', {1, 2, 2, 1, 1, []}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % t2 is s2 modified to have similar indexing to s1
    % 6
    % |
    % 5 — 3
    % |
    % 4 — 2
    % |
    % 1
    testCase.TestData.t2s1 = struct(...
        'Name', {'1', '2', '3', '', '', 'Adam'}, ...
        'parent', {4, 4, 5, 5, 6, [],}, ...
        'child', {[], [], [], [1, 2], [4, 3], 5}, ...
        'type', {0, 0, 0, 1, 2, 3}, ...
        'time', {0, 0, 0, 1, 2, 3}, ...
        'sibling', {1, 2, 2, 1, 1, []}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % s3
    % 6
    % |
    % 4 — 5 — 3
    % |   |
    % 1   2
    testCase.TestData.s3 = struct(...
        'Name', {'1', '2', '3', '', '', 'Adam'}, ...
        'parent', {4, 5, 5, 6, 4, [],}, ...
        'child', {[], [], [], [1, 5], [2, 3], 4}, ...
        'type', {0, 0, 0, 2, 1, 3}, ...
        'time', {0, 0, 0, 2, 1, 3}, ...
        'sibling', {1, 1, 2, 1, 2, []}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % t3 is s3 modified to have similar indexing to s1
    % 6
    % |
    % 5 — 4 — 3
    % |   |
    % 1   2
    testCase.TestData.t3s1 = struct(...
        'Name', {'1', '2', '3', '', '', 'Adam'}, ...
        'parent', {5, 4, 4, 5, 6, [],}, ...
        'child', {[], [], [], [2, 3], [1, 4], 5}, ...
        'type', {0, 0, 0, 1, 2, 3}, ...
        'time', {0, 0, 0, 1, 2, 3}, ...
        'sibling', {1, 1, 2, 2, 1, []}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % s4 leaf names no longer match index (leaf index.name)
    % 1
    % |
    % 2 — 4 — 5.2
    % |   |
    % 3.1 6.3
    testCase.TestData.s4 = struct(...
        'Name', {'Adam', '', '1', '', '2', '3'}, ...
        'parent', {[], 1, 2, 2, 4, 4,}, ...
        'child', {2, [3, 4], [], [6, 5], [], []}, ...
        'type', {3, 2, 0, 1, 0, 0}, ...
        'time', {3, 2, 0, 1, 0, 0}, ...
        'sibling', {[], 1, 1, 2, 2, 1}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % t4 is s4 modified to have similar indexing to s1, leaf names match indices
    % but sibling order is different from t3s1
    % 6
    % |
    % 5 — 4 — 2
    % |   |
    % 1   3
    testCase.TestData.t4s1 = struct(...
        'Name', {'1', '2', '3', '', '', 'Adam'}, ...
        'parent', {5, 4, 4, 5, 6, [],}, ...
        'child', {[], [], [], [3, 2], [1, 4], 5}, ...
        'type', {0, 0, 0, 1, 2, 3}, ...
        'time', {0, 0, 0, 1, 2, 3}, ...
        'sibling', {1, 2, 1, 2, 1, []}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % s5 is a permutation of the indices in s1, different leaf times
    % 1
    % |
    % 2 — 4.3
    % |
    % 3 — 5.1
    % |
    % 6.2
    testCase.TestData.s5 = struct(...
        'Name', {'Adam', '', '', '3', '1', '2'}, ...
        'parent', {[], 1, 2, 2, 3, 3}, ...
        'child', {2, [3, 4], [6, 5], [], [], []}, ...
        'type', {3, 2, 1, 0, 0, 0}, ...
        'time', {3, 2, 1, 0.1, 0.3, 0.2}, ...
        'sibling', {[], 1, 1, 2, 2, 1}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % t5 is similar to s1 but with different sibling orders and leaf times
    % 6
    % |
    % 5 — 3.3
    % |
    % 4 — 1.1
    % |
    % 2.2
    testCase.TestData.t5s1 = struct(...
        'Name', {'1', '2', '3', '', '', 'Adam'}, ...
        'parent', {4, 4, 5, 5, 6, []}, ...
        'child', {[], [], [], [2, 1], [4, 3], 5}, ...
        'type', {0, 0, 0, 1, 2, 3}, ...
        'time', {0.3, 0.2, 0.1, 1, 2, 3}, ...
        'sibling', {2, 1, 2, 1, 1, []}, ...
        'mark', {0, 0, 0, 0, 0, 0});

    % s6 is the true tree for 10-leaf synthetic data set with 3 clades (SIM-N)
    % s7 and s8 are the 5e4th and 1e5th iterations of a MCMC respectively
    % initialised at s6

    % % s6 ((1, (2, (4, 5))), (3, ((6, 7), (8, (9, 10))))
    % % s7 ((2, (4, 5)), (1, (3, ((6, 7), (8, (9, 10))))))
    % % s8 (1, ((2, (4, 5)), (3, ((6, 7), (8, (9, 10))))))
    %
    % (
    %     (1, (2, (4, 5))),
    %     (3, ((6, 7), (8, (9, 10))))
    % )
    %
    % ((1, (2, (4, 5))), (3, ((6, 7), (8, (9, 10))))
    % ((1, (2, (4, 5))), (3, ((6, 7), (8, (9, 10)))))
    %
    % % fStr = 'borrowing/housekeepingTestDataState%i.mat';
    % % getTree = @(i) getfield(load(sprintf(fStr, i)), 'state', 'tree');
    % % testCase.TestData.s6 = getTree(6);
    % % testCase.TestData.s7 = getTree(7);
    % % testCase.TestData.s8 = getTree(8);

    % 1 an outgroup of left subtree
    testCase.TestData.s6 = rnextree(...
        '((1, (2, (4, 5))), (3, ((6, 7), (8, (9, 10)))))');
    % 1 an outgroup of right subtree
    testCase.TestData.s7 = rnextree(...
        '((2, (4, 5)), (1, (3, ((6, 7), (8, (9, 10))))))');
    % 1 an overall outgroup
    testCase.TestData.s8 = rnextree(...
        '(1, ((2, (4, 5)), (3, ((6, 7), (8, (9, 10))))))');

    % Indices of s7 switched to same roles in s6
    testCase.TestData.t7s6 = struct(...
        'Name', {[], [], '1', [], '2', [], '4', '5', [], '3', [], [], '6', '7', [], '8', [], '9', '10', 'Adam'}, ...
        'parent', {20, 1, 2, 1, 4, 4, 6, 6, 2, 9, 9, 11, 12, 12, 11, 15, 15, 17, 17, []}, ...
        'child', {[4, 2], [3, 9], [], [5, 6], [], [7, 8], [], [], [10, 11], [], [12, 15], [13, 14], [], [], [16, 17], [], [18, 19], [], [], 1}, ...
        'type', {2, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 3}, ...
        'time', {6, 5, 4, 5, 4, 4, 3, 3, 4, 3, 3, 2, 1, 1, 2, 1, 1, 0, 0, testCase.TestData.s7(20).time}, ...
        'sibling', {1, 2, 1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, []}, ...
        'mark', {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    % Much easier problem as only pa(1) and MRCA(2, ..., 10) need to switch
    testCase.TestData.t8s6 = struct(...
        'Name', {[], [], '1', [], '2', [], '4', '5', [], '3', [], [], '6', '7', [], '8', [], '9', '10', 'Adam'}, ...
        'parent', {20, 1, 1, 2, 4, 4, 6, 6, 2, 9, 9, 11, 12, 12, 11, 15, 15, 17, 17, []}, ...
        'child', {[3, 2], [4, 9], [], [5, 6], [], [7, 8], [], [], [10, 11], [], [12, 15], [13, 14], [], [], [16, 17], [], [18, 19], [], [], 1}, ...
        'type', {2, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 3}, ...
        'time', {6, 5, 5, 4, 3, 3, 2, 2, 4, 3, 3, 2, 1, 1, 2, 1, 1, 0, 0, testCase.TestData.s7(20).time}, ...
        'sibling', {1, 2, 1, 1, 1, 2, 1, 2, 2,  1, 2, 1, 1, 2, 2, 1, 2, 1, 2, []}, ...
        'mark', {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
end

function teardownOnce(testCase)
    rmpath('debugging', '.');
end

function tests = housekeepingTest
    % Unit-testing likelihood calculations in patternCounts
    tests = functiontests(localfunctions);
end

%% Test Functions
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


% Helper functions
function checkHousekeeping(testCase, tObs, tExp)
    fn = {'Name', 'parent', 'child', 'type', 'time'}'; % fieldnames(tExp);
    for j = find(cellfun(@(x) ~isequal(x, 'time'), fn))'
        assertTrue(testCase, ...
                   all(cellfun(@(a, b) isequaln(sort(a), sort(b)), ...
                               {tObs.(fn{j})}, ...
                               {tExp.(fn{j})})));
    end
end

% function checkHousekeeping(testCase, s1, t1)
%     % Check all fields which should match do (we ignore time)
%     fn = fieldnames(s1);
%     t2 = housekeeping(s1, t1);
%     [c1, ~] = ismember(getLeafArray(s1), getLeafArray(t2), 'rows');
%     for j = find(cellfun(@(x) ~isequal(x, 'time'), fn))'
%         fs1 = {s1.(fn{j})};
%         ft2 = {t2.(fn{j})};
%         assertTrue(testCase, all(cellfun(@isequaln, fs1(~c1), ft2(~c1))));
%     end
% end

% Optional file fixtures
function setupOnce(testCase)
    % Global variables
    GlobalSwitches;
    GlobalValues;

    % Fake trees
    % Root and left child are sibling 1, right child is sibling 2
    % Leaves are at time 0, internal node at 1, root at 2 and Adam at 3

    % Before housekeeping

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
        'sibling', {1, 2, 2, 1, 1, []});

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
        'sibling', {1, 2, 2, 1, 1, []});

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
        'sibling', {1, 1, 2, 1, 2, []});

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
        'sibling', {[], 1, 1, 2, 2, 1});

    % After housekeeping

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
        'sibling', {1, 2, 2, 1, 1, []});

    % t3 is s3 modified to have similar indexing to
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
        'sibling', {1, 1, 2, 2, 1, []});

    % t4 is s4 modified to have similar indexing to s1, leaf names match indices
    % 6
    % |
    % 5 — 4 — 3
    % |   |
    % 1   2
    testCase.TestData.t4s1 = testCase.TestData.t3s1;
end

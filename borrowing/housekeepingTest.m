function tests = housekeepingTest
    % Unit-testing likelihood calculations in patternCounts
    tests = functiontests(localfunctions);
end

%% Test Functions
function s1s1Test(testCase)
    % Identical inputs yield identical outputs, no catastrophes
    prior = testCase.TestData.p1;
    state1 = partialMakeState(prior, testCase.TestData.s1);
    nstate1 = housekeeping(state1, state1);
    checkHousekeeping(testCase, nstate1, state1);
end

function s1s2Test(testCase)
    % Switch root and internal node, overlapping subtree sets, catastrophe on s2
    prior = testCase.TestData.p1;
    state1 = partialMakeState(prior, testCase.TestData.s1);
    s2 = testCase.TestData.s2;
    s2(5).catloc = 0.5;
    state2 = partialMakeState(prior, s2);
    nstate2Obs = housekeeping(state1, state2);
    t2s1 = testCase.TestData.t2s1;
    t2s1(4).catloc = 0.5;
    nstate2Exp = partialMakeState(prior, t2s1);
    checkHousekeeping(testCase, nstate2Obs, nstate2Exp);
end

function s1s3Test(testCase)
    % Switch root and internal node, non-overlapping subtree sets, no clades or
    % catastrophes
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state3 = partialMakeState(prior, testCase.TestData.s3);
    nstate3Obs = housekeeping(state1, state3);
    nstate3Exp = partialMakeState(prior, testCase.TestData.t3s1);
    checkHousekeeping(testCase, nstate3Obs, nstate3Exp);
end

function s1s4Test(testCase)
    % Switch Adam, root internal node and leaves, non-overlapping subtree sets,
    % no clades or catastrophes
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state4 = partialMakeState(prior, testCase.TestData.s4);
    nstate4Obs = housekeeping(state1, state4);
    nstate4Exp = partialMakeState(prior, testCase.TestData.t4s1);
    checkHousekeeping(testCase, nstate4Obs, nstate4Exp);
end

function s1s5Test(testCase)
    % Permutation of indices, leaves at different times, no clades or
    % catastrophes
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state5 = partialMakeState(prior, testCase.TestData.s5);
    nstate5Obs = housekeeping(state1, state5);
    nstate5Exp = partialMakeState(prior, testCase.TestData.t5s1);
    checkHousekeeping(testCase, nstate5Obs, nstate5Exp);
end

function s6s7Test(testCase)
    % Permutation of indices, leaves at different times, no clades or
    % catastrophes
    prior = testCase.TestData.p6;
    state6 = partialMakeState(prior, testCase.TestData.s6);
    state7 = partialMakeState(prior, testCase.TestData.s7);
    nstate7Obs = housekeeping(state6, state7);
    nstate7Exp = partialMakeState(prior, testCase.TestData.t7s6);
    checkHousekeeping(testCase, nstate7Obs, nstate7Exp);
end

function s6s8Test(testCase)
    % Only one node swap, leaves at different times, no clades or catastrophes
    prior = testCase.TestData.p6;
    state6 = partialMakeState(prior, testCase.TestData.s6);
    state8 = partialMakeState(prior, testCase.TestData.s8);
    nstate8Obs = housekeeping(state6, state8);
    nstate8Exp = partialMakeState(prior, testCase.TestData.t8s6);
    checkHousekeeping(testCase, nstate8Obs, nstate8Exp);
end

% Helper functions
function checkHousekeeping(testCase, stateObs, stateExp)
    % Compare tree variables
    tfn = {'Name', 'parent', 'child', 'type', 'time', 'sibling'};
    sObs = stateObs.tree;
    sExp = stateExp.tree;
    for j = 1:length(tfn)
        assertTrue(testCase, all(cellfun(@(a, b) isequaln(a, b), ...
                                         {sObs.(tfn{j})}, {sExp.(tfn{j})})));
    end
    assertTrue(testCase, equaltrees(sObs, sExp) == 1);
    % Compare state variables
    sfn = {'leaves', 'root', 'nodes', 'cat', 'claderoot'};
    for j = 1:length(sfn)
        assertTrue(testCase, all(stateObs.(sfn{j}) == stateExp.(sfn{j})));
    end
end

function state = partialMakeState(prior, s)
    % Same as makestate but without data and associated calculations
    % Assume catastrophes are already on tree
    state = tree2state(s);

    state.claderoot = [];
    if prior.isclade
       state.tree = treeclades(state.tree, prior.clade);
       state = UpdateClades(state, [state.leaves, state.nodes], ...
                            size(prior.clade, 2));
    end

    state.cat = cellfun('length', {state.tree.catloc});
    state.cat = state.cat(:);
    state.ncat = sum(state.cat);
    state.length = TreeLength(state.tree, state.root);
end

function prior = popClades()
    prior = getfield(pop('model'), 'prior');
end

% Optional file fixtures
function setupOnce(testCase)
    % Setting up
    GlobalSwitches;
    GlobalValues;

    % Root and left child are sibling 1, right child is sibling 2
    % Leaves are at time 0, internal node at 1, root at 2 and Adam at 3
    % TODO: Most of these trees are only used once so move inside test functions

    emptyTreeStruct = @(n) repmat(TreeNode([], [], [], [], [], []), 1, n);

    % s1
    % 6
    % |
    % 5 — 3
    % |
    % 4 — 2
    % |
    % 1
    s1 = emptyTreeStruct(6);
    [s1.Name] = deal('1', '2', '3', '', '', 'Adam');
    [s1.parent] = deal(4, 4, 5, 5, 6, []);
    [s1.child] = deal([], [], [], [1, 2], [4, 3], 5);
    [s1.type] = deal(0, 0, 0, 1, 2, 3);
    [s1.time] = deal(0, 0, 0, 1, 2, 3);
    [s1.sibling] = deal(1, 2, 2, 1, 1, []);
    testCase.TestData.s1 = s1;

    % s2
    % 6
    % |
    % 4 — 3
    % |
    % 5 — 2
    % |
    % 1
    s2 = emptyTreeStruct(6);
    [s2.Name] = deal('1', '2', '3', '', '', 'Adam');
    [s2.parent] = deal(5, 5, 4, 6, 4, []);
    [s2.child] = deal([], [], [], [5, 3], [1, 2], 4);
    [s2.type] = deal(0, 0, 0, 2, 1, 3);
    [s2.time] = deal(0, 0, 0, 2, 1, 3);
    [s2.sibling] = deal(1, 2, 2, 1, 1, []);
    testCase.TestData.s2 = s2;

    % t2s1 is s2 modified to have similar indexing to s1
    % 6
    % |
    % 5 — 3
    % |
    % 4 — 2
    % |
    % 1
    t2s1 = emptyTreeStruct(6);
    [t2s1.Name] = deal('1', '2', '3', '', '', 'Adam');
    [t2s1.parent] = deal(4, 4, 5, 5, 6, []);
    [t2s1.child] = deal([], [], [], [1, 2], [4, 3], 5);
    [t2s1.type] = deal(0, 0, 0, 1, 2, 3);
    [t2s1.time] = deal(0, 0, 0, 1, 2, 3);
    [t2s1.sibling] = deal(1, 2, 2, 1, 1, []);
    testCase.TestData.t2s1 = t2s1;

    % s3
    % 6
    % |
    % 4 — 5 — 3
    % |   |
    % 1   2
    s3 = emptyTreeStruct(6);
    [s3.Name] = deal('1', '2', '3', '', '', 'Adam');
    [s3.parent] = deal(4, 5, 5, 6, 4, []);
    [s3.child] = deal([], [], [], [1, 5], [2, 3], 4);
    [s3.type] = deal(0, 0, 0, 2, 1, 3);
    [s3.time] = deal(0, 0, 0, 2, 1, 3);
    [s3.sibling] = deal(1, 1, 2, 1, 2, []);
    testCase.TestData.s3 = s3;

    % t3s1 is s3 modified to have similar indexing to s1
    % 6
    % |
    % 5 — 4 — 3
    % |   |
    % 1   2
    t3s1 = emptyTreeStruct(6);
    [t3s1.Name] = deal('1', '2', '3', '', '', 'Adam');
    [t3s1.parent] = deal(5, 4, 4, 5, 6, []);
    [t3s1.child] = deal([], [], [], [2, 3], [1, 4], 5);
    [t3s1.type] = deal(0, 0, 0, 1, 2, 3);
    [t3s1.time] = deal(0, 0, 0, 1, 2, 3);
    [t3s1.sibling] = deal(1, 1, 2, 2, 1, []);
    testCase.TestData.t3s1 = t3s1;

    % s4 leaf names no longer match index (leaf index.name)
    % 1
    % |
    % 2 — 4 — 5.2
    % |   |
    % 3.1 6.3
    s4 = emptyTreeStruct(6);
    [s4.Name] = deal('Adam', '', '1', '', '2', '3');
    [s4.parent] = deal([], 1, 2, 2, 4, 4);
    [s4.child] = deal(2, [3, 4], [], [6, 5], [], []);
    [s4.type] = deal(3, 2, 0, 1, 0, 0);
    [s4.time] = deal(3, 2, 0, 1, 0, 0);
    [s4.sibling] = deal([], 1, 1, 2, 2, 1);
    testCase.TestData.s4 = s4;

    % t4s1 is s4 modified to have similar indexing to s1, leaf names match
    % indices but sibling order is different from t3s1
    % 6
    % |
    % 5 — 4 — 2
    % |   |
    % 1   3
    t4s1 = emptyTreeStruct(6);
    [t4s1.Name] = deal('1', '2', '3', '', '', 'Adam');
    [t4s1.parent] = deal(5, 4, 4, 5, 6, []);
    [t4s1.child] = deal([], [], [], [3, 2], [1, 4], 5);
    [t4s1.type] = deal(0, 0, 0, 1, 2, 3);
    [t4s1.time] = deal(0, 0, 0, 1, 2, 3);
    [t4s1.sibling] = deal(1, 2, 1, 2, 1, []);
    testCase.TestData.t4s1 = t4s1;

    % s5 is a permutation of the indices in s1, different leaf times
    % 1
    % |
    % 2 — 4.3
    % |
    % 3 — 5.1
    % |
    % 6.2
    s5 = emptyTreeStruct(6);
    [s5.Name] = deal('Adam', '', '', '3', '1', '2');
    [s5.parent] = deal([], 1, 2, 2, 3, 3);
    [s5.child] = deal(2, [3, 4], [6, 5], [], [], []);
    [s5.type] = deal(3, 2, 1, 0, 0, 0);
    [s5.time] = deal(3, 2, 1, 0.1, 0.3, 0.2);
    [s5.sibling] = deal([], 1, 1, 2, 2, 1);
    testCase.TestData.s5 = s5;

    % t5s1 is similar to s1 but with different sibling orders and leaf times
    % 6
    % |
    % 5 — 3.3
    % |
    % 4 — 1.1
    % |
    % 2.2
    t5s1 = emptyTreeStruct(6);
    [t5s1.Name] = deal('1', '2', '3', '', '', 'Adam');
    [t5s1.parent] = deal(4, 4, 5, 5, 6, []);
    [t5s1.child] = deal([], [], [], [2, 1], [4, 3], 5);
    [t5s1.type] = deal(0, 0, 0, 1, 2, 3);
    [t5s1.time] = deal(0.3, 0.2, 0.1, 1, 2, 3);
    [t5s1.sibling] = deal(2, 1, 2, 1, 1, []);
    testCase.TestData.t5s1 = t5s1;

    % Clade constraint for previous trees
    p1 = popClades();
    p1.isclade = 1;
    p1.clade = {pop('clade')};
    p1.clade{1}.language = {'1', '2', '3'};
    testCase.TestData.p1 = p1;

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
    t7s6 = emptyTreeStruct(20);
    [t7s6.Name] = deal([], [], '1', [], '2', [], '4', '5', [], '3', [], [], '6', '7', [], '8', [], '9', '10', 'Adam');
    [t7s6.parent] = deal(20, 1, 2, 1, 4, 4, 6, 6, 2, 9, 9, 11, 12, 12, 11, 15, 15, 17, 17, []);
    [t7s6.child] = deal([4, 2], [3, 9], [], [5, 6], [], [7, 8], [], [], [10, 11], [], [12, 15], [13, 14], [], [], [16, 17], [], [18, 19], [], [], 1);
    [t7s6.type] = deal(2, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 3);
    [t7s6.time] = deal(6, 5, 4, 5, 4, 4, 3, 3, 4, 3, 3, 2, 1, 1, 2, 1, 1, 0, 0, testCase.TestData.s7(20).time);
    [t7s6.sibling] = deal(1, 2, 1, 1, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, []);
    testCase.TestData.t7s6 = t7s6;

    % Much easier problem as only pa(1) and MRCA(2, ..., 10) need to switch
    t8s6 = emptyTreeStruct(20);
    [t8s6.Name] = deal([], [], '1', [], '2', [], '4', '5', [], '3', [], [], '6', '7', [], '8', [], '9', '10', 'Adam');
    [t8s6.parent] = deal(20, 1, 1, 2, 4, 4, 6, 6, 2, 9, 9, 11, 12, 12, 11, 15, 15, 17, 17, []);
    [t8s6.child] = deal([3, 2], [4, 9], [], [5, 6], [], [7, 8], [], [], [10, 11], [], [12, 15], [13, 14], [], [], [16, 17], [], [18, 19], [], [], 1);
    [t8s6.type] = deal(2, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 3);
    [t8s6.time] = deal(6, 5, 5, 4, 3, 3, 2, 2, 4, 3, 3, 2, 1, 1, 2, 1, 1, 0, 0, testCase.TestData.s7(20).time);
    [t8s6.sibling] = deal(1, 2, 1, 1, 1, 2, 1, 2, 2,  1, 2, 1, 1, 2, 2, 1, 2, 1, 2, []);
    testCase.TestData.t8s6 = t8s6;

    % Clade constraints for previous trees
    % Clades 2 and 3 aren't used for any moves currently
    p6 = popClades();
    p6.isclade = 1;
    p6.clade = repmat({pop('clade')}, 1, 3);
    p6.clade{1}.language = {'1'};
    p6.clade{2}.language = {'8'};
    p6.clade{3}.language = {'6', '7', '8', '9', '10'};
    testCase.TestData.p6 = p6;
end

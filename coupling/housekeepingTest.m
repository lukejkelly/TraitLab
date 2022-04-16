function tests = housekeepingTest
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

function s1s9Test(testCase)
    % Only root and adam node indices differ
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state9 = partialMakeState(prior, testCase.TestData.s9);
    nstate9Obs = housekeeping(state1, state9);
    nstate9Exp = partialMakeState(prior, testCase.TestData.t9s1);
    checkHousekeeping(testCase, nstate9Obs, nstate9Exp);
end

function s1s10Test(testCase)
    % Only adam and internal node indices differ
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state10 = partialMakeState(prior, testCase.TestData.s10);
    nstate10Obs = housekeeping(state1, state10);
    nstate10Exp = partialMakeState(prior, testCase.TestData.t10s1);
    checkHousekeeping(testCase, nstate10Obs, nstate10Exp);
end

function s1s11Test(testCase)
    % Only adam and leaf node indices differ
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state11 = partialMakeState(prior, testCase.TestData.s11);
    nstate11Obs = housekeeping(state1, state11);
    nstate11Exp = partialMakeState(prior, testCase.TestData.t11s1);
    checkHousekeeping(testCase, nstate11Obs, nstate11Exp);
end

function s1s12Test(testCase)
    % Only root and internal node indices differ
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state12 = partialMakeState(prior, testCase.TestData.s12);
    nstate12Obs = housekeeping(state1, state12);
    nstate12Exp = partialMakeState(prior, testCase.TestData.t12s1);
    checkHousekeeping(testCase, nstate12Obs, nstate12Exp);
end

function s1s13Test(testCase)
    % Only root and leaf node indices differ
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s1);
    state13 = partialMakeState(prior, testCase.TestData.s13);
    nstate13Obs = housekeeping(state1, state13);
    nstate13Exp = partialMakeState(prior, testCase.TestData.t13s1);
    checkHousekeeping(testCase, nstate13Obs, nstate13Exp);
end

function s14t14Test(testCase)
    % Two common clades, check labelling within
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s14_1);
    state2 = partialMakeState(prior, testCase.TestData.s14_2);
    nstate2Obs = housekeeping(state1, state2);
    nstate2Exp = partialMakeState(prior, testCase.TestData.t14);
    checkHousekeeping(testCase, nstate2Obs, nstate2Exp);
end

function s15t15Test(testCase)
    % Two common clades, check labelling within
    prior = popClades();
    state1 = partialMakeState(prior, testCase.TestData.s15_1);
    state2 = partialMakeState(prior, testCase.TestData.s15_2);
    nstate2Obs = housekeeping(state1, state2);
    nstate2Exp = partialMakeState(prior, testCase.TestData.t15);
    checkHousekeeping(testCase, nstate2Obs, nstate2Exp);
end

% Helper functions
function checkHousekeeping(testCase, stateObs, stateExp)
    % Compare tree variables
    sObs = stateObs.tree;
    sExp = stateExp.tree;
    for tfn = {'Name', 'parent', 'child', 'type', 'time', 'sibling'}
        assertEqual(testCase, {sObs.(tfn{:})}, {sExp.(tfn{:})});
    end
    assertEqual(testCase, equaltrees(sObs, sExp), 1);
    % Compare state variables
    sfn = {'leaves', 'root', 'nodes', 'cat', 'claderoot'};
    for j = 1:length(sfn)
        assertEqual(testCase, stateObs.(sfn{j}), stateExp.(sfn{j}));
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

    checktree(state.tree, state.NS);
end

function prior = popClades()
    prior = getfield(pop('model'), 'prior');
end

% Optional file fixtures
function setupOnce(testCase)
    % Setting up
    GlobalSwitches;
    GlobalValues;

    global BORROWING
    testCase.TestData.BORROWING = BORROWING;
    BORROWING = 0;

    emptyTreeStruct = @(n) repmat(TreeNode([], [], [], [], [], []), 1, n);

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

    % t5s1 is similar to s1 but with different leaf times
    % 6
    % |
    % 5 — 3.3
    % |
    % 4 — 2.2
    % |
    % 1.1
    t5s1 = emptyTreeStruct(6);
    [t5s1.Name] = deal('1', '2', '3', '', '', 'Adam');
    [t5s1.parent] = deal(4, 4, 5, 5, 6, []);
    [t5s1.child] = deal([], [], [], [1, 2], [4, 3], 5);
    [t5s1.type] = deal(0, 0, 0, 1, 2, 3);
    [t5s1.time] = deal(0.3, 0.2, 0.1, 1, 2, 3);
    [t5s1.sibling] = deal(1, 2, 2, 1, 1, []);
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

    % Only Adam and root differ
    % Cannot plot this tree as draw assumes Adam is last index in struct

    % s9
    % 5
    % |
    % 6 — 3
    % |
    % 4 — 2
    % |
    % 1
    s9 = emptyTreeStruct(6);
    [s9.Name] = deal('1', '2', '3', '', 'Adam', '');
    [s9.parent] = deal(4, 4, 6, 6, [], 5);
    [s9.child] = deal([], [], [], [1, 2], 6, [4, 3]);
    [s9.type] = deal(0, 0, 0, 1, 3, 2);
    [s9.time] = deal(0, 0, 0, 1, 3, 2);
    [s9.sibling] = deal(1, 2, 2, 1, [], 1);
    testCase.TestData.s9 = s9;

    % t9s1 is s1 exactly
    testCase.TestData.t9s1 = s1;

    % Only adam and internal node differ
    % s10
    % 4
    % |
    % 5 — 3
    % |
    % 6 — 2
    % |
    % 1
    s10 = emptyTreeStruct(6);
    [s10.Name] = deal('1', '2', '3', 'Adam', '', '');
    [s10.parent] = deal(6, 6, 5, [], 4, 5);
    [s10.child] = deal([], [], [], 5, [6, 3], [1, 2]);
    [s10.type] = deal(0, 0, 0, 3, 2, 1);
    [s10.time] = deal(0, 0, 0, 3, 2, 1);
    [s10.sibling] = deal(1, 2, 2, [], 1, 1);
    testCase.TestData.s10 = s10;

    % t10s1 is s1 exactly
    testCase.TestData.t10s1 = s1;

    % Only adam and leaf differ
    % s11
    % 2
    % |
    % 5 — 3
    % |
    % 4 — 6.2
    % |
    % 1
    s11 = emptyTreeStruct(6);
    [s11.Name] = deal('1', 'Adam', '3', '', '', '2');
    [s11.parent] = deal(4, [], 5, 5, 2, 4);
    [s11.child] = deal([], 5, [], [1, 6], [4, 3], []);
    [s11.type] = deal(0, 3, 0, 1, 2, 0);
    [s11.time] = deal(0, 3, 0, 1, 2, 0);
    [s11.sibling] = deal(1, [], 2, 1, 1, 2);
    testCase.TestData.s11 = s11;

    % t11s1 is s1 exactly
    testCase.TestData.t11s1 = s1;

    % Only root and internal node differ
    % s12
    % 6
    % |
    % 4 — 3
    % |
    % 5 — 2
    % |
    % 1
    s12 = emptyTreeStruct(6);
    [s12.Name] = deal('1', '2', '3', '', '', 'Adam');
    [s12.parent] = deal(5, 5, 4, 6, 4, []);
    [s12.child] = deal([], [], [], [5, 3], [1, 2], 4);
    [s12.type] = deal(0, 0, 0, 2, 1, 3);
    [s12.time] = deal(0, 0, 0, 2, 1, 3);
    [s12.sibling] = deal(1, 2, 2, 1, 1, []);
    testCase.TestData.s12 = s12;

    % t11s1 is s1 exactly
    testCase.TestData.t12s1 = s1;

    % Only root and leaf differ
    % s13
    % 6
    % |
    % 2 — 3
    % |
    % 4 — 5.2
    % |
    % 1
    s13 = emptyTreeStruct(6);
    [s13.Name] = deal('1', '', '3', '', '2', 'Adam');
    [s13.parent] = deal(4, 6, 2, 2, 4, []);
    [s13.child] = deal([], [4, 3], [], [1, 5], [], 2);
    [s13.type] = deal(0, 2, 0, 1, 0, 3);
    [s13.time] = deal(0, 2, 0, 1, 0, 3);
    [s13.sibling] = deal(1, 1, 2, 1, 2, []);
    testCase.TestData.s13 = s13;

    % t11s1 is s1 exactly
    testCase.TestData.t13s1 = s1;

    % Checking within-clade labelling
    % Afterwards, left and right subtrees in s2 should use same index sets as s1

    %              16
    %               |
    %         ______9_________
    %    ___10____       ___13___
    % __11_   __12_   __14_   __15_
    % 1   2   3   4   5   6   7   8
    s14_1 = emptyTreeStruct(16);
    [s14_1.Name] = deal('1', '2', '3', '4', '5', '6', '7', '8', '', '', '', '', '', '', '', 'Adam');
    [s14_1.parent] = deal(11, 11, 12, 12, 14, 14, 15, 15, 16, 9, 10, 10, 9, 13, 13, []);
    [s14_1.child] = deal([], [], [], [], [], [], [], [], [10, 13], [11, 12], [1, 2], [3, 4], [14, 15], [5, 6], [7, 8], 9);
    [s14_1.type] = deal(0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 3);
    [s14_1.time] = deal(0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 2, 1, 1, 4);
    [s14_1.sibling] = deal(1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2, 1, 2, []);
    testCase.TestData.s14_1 = s14_1;

    %              16
    %               |
    %         _____11__________
    %     ___9____        ___12____
    % ___15___    |   ___10___    |
    % |   __13_   |   |   __14_   |
    % 1   2   3   4   5   6   7   8
    s14_2 = emptyTreeStruct(16);
    [s14_2.Name] = deal('1', '2', '3', '4', '5', '6', '7', '8', '', '', '', '', '', '', '', 'Adam');
    [s14_2.parent] = deal(15, 13, 13, 9, 10, 14, 14, 12, 11, 12, 16, 11, 15, 10, 9, []);
    [s14_2.child] = deal([], [], [], [], [], [], [], [], [15, 4], [5, 14], [9, 12], [10, 8], [2, 3], [6, 7], [1, 13], 11);
    [s14_2.type] = deal(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 1, 1, 3);
    [s14_2.time] = deal(0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 4, 3, 1, 1, 2, 5);
    [s14_2.sibling] = deal(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 1, []);
    testCase.TestData.s14_2 = s14_2;

    %              16
    %               |
    %         ______9__________
    %     ___10____        ___13___
    % ___11___    |   ___15___    |
    % |   __12_   |   |   __14_   |
    % 1   2   3   4   5   6   7   8
    t14 = emptyTreeStruct(16);
    [t14.Name] = deal('1', '2', '3', '4', '5', '6', '7', '8', '', '', '', '', '', '', '', 'Adam');
    [t14.parent] = deal(11, 12, 12, 10, 15,  14, 14, 13, 16, 9, 10, 11, 9, 15, 13, []);
    [t14.child] = deal([], [], [], [], [], [], [], [], [10, 13], [11, 4], [1, 12], [2, 3], [15, 8], [6, 7], [5, 14], 9);
    [t14.type] = deal(0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 1, 1, 1, 3);
    [t14.time] = deal(0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 2, 1, 3, 1, 2, 5);
    [t14.sibling] = deal(1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 2, 2, 2, 1, []);
    testCase.TestData.t14 = t14;

    % Same as 14 but on a bigger scale
    s15_1 = emptyTreeStruct(24);
    [s15_1.Name] = deal('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '', '', '', '', '', '', '', '', '', '', '', 'Adam');
    [s15_1.parent] = deal(23, 23, 22, 21, 20, 19, 17, 16, 15, 14, 13, 13, 14, 15, 16, 17, 18, 24, 18, 19, 20, 21, 22, []);
    [s15_1.child] = deal([], [], [], [], [], [], [], [], [], [], [], [], [11, 12], [10, 13], [9, 14], [8, 15], [7, 16], [19, 17], [20, 6], [21, 5], [22, 4], [23, 3], [1, 2], 18);
    [s15_1.type] = deal(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 3);
    [s15_1.time] = deal(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1, 7);
    [s15_1.sibling] = deal(1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, []);
    testCase.TestData.s15_1 = s15_1;

    s15_2 = emptyTreeStruct(24);
    [s15_2.Name] = deal('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '', '', '', '', '', '', '', '', '', '', '', 'Adam');
    [s15_2.parent] = deal(17, 16, 15, 14, 13, 13, 18, 18, 19, 20, 21, 22, 14, 15, 16, 17, 23, 19, 20, 21, 22, 23, 24, []);
    [s15_2.child] = deal([], [], [], [], [], [], [], [], [], [], [], [], [5, 6], [4, 13], [3, 14], [2, 15], [1, 16], [7, 8], [18, 9], [19, 10], [20, 11], [21, 12], [17, 22], 23);
    [s15_2.type] = deal(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3);
    [s15_2.time] = deal(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4,  5, 1, 2, 3, 4, 5, 6, 7);
    [s15_2.sibling] = deal(1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, []);
    testCase.TestData.s15_2 = s15_2;

    t15 = emptyTreeStruct(24);
    [t15.Name] = deal('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '', '', '', '', '', '', '', '', '', '', '', 'Adam');
    [t15.parent] = deal(19, 23, 22, 21, 20, 20, 16,  16,  15, 13, 14, 17, 14, 17, 13, 15, 18, 24,  18, 21, 22, 23, 19, []);
    [t15.child] = deal([], [], [], [], [], [], [], [], [], [], [], [], [15, 10], [13, 11], [16, 9], [7, 8], [14,  12], [19, 17], [1, 23], [5, 6], [4, 20], [3, 21], [2, 22], 18);
    [t15.type] = deal(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 3);
    [t15.time] = deal(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 4, 2, 1, 5, 6, 5, 1, 2, 3, 4, 7);
    [t15.sibling] = deal(1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 2, []);
    testCase.TestData.t15 = t15;
end

function teardownOnce(testCase)
    global BORROWING
    BORROWING = testCase.TestData.BORROWING;
end

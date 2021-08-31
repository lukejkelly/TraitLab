function tests = GetLegalCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    rangeL = 6:10;
    n_i = length(rangeL);
    n_j = 1e2;
    [new_x, new_y, q1_x, q1_y, q2_x, q2_y] = deal(zeros(n_i, n_j));
    for i = 1:n_i
        [state_x, state_y] = coupledStates(rangeL(i));
        s_x = state_x.tree;
        s_y = state_y.tree;
        root = state_x.root;
        for j = 1:n_j
            [old_x, old_y] = sampleBranchByCatCountCoupled(state_x, state_y);
            assertEqual(testCase, old_x, old_y);
            [new_x(i, j), new_y(i, j), q1_x(i, j), q1_y(i, j), q2_x(i, j), ...
                q2_y(i, j)] = GetLegalCoupled(s_x, s_y, old_x, old_y, root);
            [new_z, q1_z, q2_z] = GetLegal(s_x, old_x, root);
            if new_x(i, j) == new_z
                assertEqual(testCase, new_x(i, j), new_z);
                assertEqual(testCase, q1_x(i, j), q1_z);
                assertEqual(testCase, q2_x(i, j), q2_z);
            end
        end
    end
    assertEqual(testCase, new_x, new_y);
    assertEqual(testCase, q1_x, q1_y);
    assertEqual(testCase, q2_x, q2_y);
end

function housekeptTest(testCase)
    rangeL = 6:10;
    n_i = length(rangeL);
    n_j = 1e4;
    [cObs, cExp] = deal(zeros(n_i, 1));
    for i = 1:n_i
        [state_x, state_y] = housekeptStates(rangeL(i));
        s_x = state_x.tree;
        s_y = state_y.tree;
        root = state_x.root;
        for j = 1:n_j
            [old_x, old_y] = sampleBranchByCatCountCoupled(state_x, state_y);
            [new_x, new_y] = GetLegalCoupled(s_x, s_y, old_x, old_y, root);
            if old_x == old_y && new_x == new_y
                cObs(i) = cObs(i) + 1;
            end
        end
        cExp(i) = getOverlap(state_x, state_y);
    end
    cObs = cObs / n_j;
    fprintf('Proportion of matching samples in each of %g trials\n', n_j);
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs - cExp));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function [state_x, state_y] = housekeptStates(L)
    [state_x, state_y] = unitTests.housekeptStates(L, 1e-2);
    for c = 1:(1 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    for c = 1:(1 + poissrnd(3))
        state_x = AddCat(state_x);
        state_y = AddCat(state_y);
    end
end

function [state_x, state_y] = coupledStates(L)
    state_x = unitTests.dummyState(ExpTree(L, 1e-2));
    for c = 1:(1 + poissrnd(3))
        state_x = AddCat(state_x);
    end
    state_y = state_x;
end

function ol = getOverlap(state_x, state_y)
    % Probability that destination branches are identical
    p_x = state_x.cat' / state_x.ncat;
    p_y = state_y.cat' / state_y.ncat;
    s_x = state_x.tree;
    s_y = state_y.tree;
    root = state_x.root;
    ol = 0;
    for i = find(p_x > 0 & p_y > 0)
        [poss_x, q_x] = GetLegalCoupled.getPoss(s_x, i, root);
        [poss_y, q_y] = GetLegalCoupled.getPoss(s_y, i, root);
        z = length(intersect(poss_x, poss_y)) / max(q_x, q_y);
        ol = ol + min(p_x(i), p_y(i)) * z;
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 0;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function tests = sampleBranchProportionalToLengthCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global ROOT
    rangeL = 6:10;
    n_i = length(rangeL);
    n_j = 1e2;
    [i_x, i_y] = deal(zeros(n_i, n_j));
    for i = 1:n_i
        state.NS = rangeL(i);
        s = ExpTree(state.NS, 1e-2);
        state.tree = s;
        state.root = find([s.type] == ROOT);
        state.length = TreeLength(s, state.root);
        for j = 1:n_j
            [i_x(i, j), i_y(i, j)] ...
                = sampleBranchProportionalToLengthCoupled(state, state);
        end
    end
    assertEqual(testCase, i_x, i_y);
end

function housekeptTest(testCase)
    rangeL = 6:10;
    n_i = length(rangeL);
    n_j = 1e4;
    [cObs, cExp] = deal(zeros(n_i, 1));
    for i = 1:n_i
        [state_x, state_y] ...
            = RscaleSubTreeCoupled.unitTests.housekeptStates(rangeL(i), 1e-2);
        for j = 1:n_j
            [i_x, i_y] = sampleBranchProportionalToLengthCoupled(state_x, ...
                                                                 state_y);
            if i_x == i_y
                cObs(i) = cObs(i) + 1;
            end
        end
        p_x = getBranchLengths(state_x) ./ state_x.length;
        p_y = getBranchLengths(state_y) ./ state_y.length;
        cExp(i) = sum(min(p_x, p_y));
    end
    cObs = cObs / n_j;
    fprintf('Proportion of matching samples in each of %g trials\n', n_j);
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs - cExp));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

% Setup and teardown functions
function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    global BORROWING;
    testCase.TestData.BORROWING = BORROWING;
    BORROWING = 1;
end

function teardownOnce(testCase)
    global BORROWING;
    BORROWING = testCase.TestData.BORROWING;
end

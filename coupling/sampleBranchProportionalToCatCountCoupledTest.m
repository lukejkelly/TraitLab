function tests = sampleBranchProportionalToCatCountCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    rangeL = 6:10;
    n_i = length(rangeL);
    n_j = 1e2;
    [i_x, i_y] = deal(zeros(n_i, n_j));
    for i = 1:n_i
        [state_x, state_y] = coupledStates(rangeL(i));
        for j = 1:n_j
            [i_x(i, j), i_y(i, j)] ...
                = sampleBranchProportionalToCatCountCoupled(state_x, state_y);
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
        [state_x, state_y] = housekeptStates(rangeL(i));
        for j = 1:n_j
            [i_x, i_y] = sampleBranchProportionalToCatCountCoupled(state_x, ...
                                                                   state_y);
            if i_x == i_y
                cObs(i) = cObs(i) + 1;
            end
        end
        p_x = state_x.cat ./ state_x.ncat;
        p_y = state_y.cat ./ state_y.ncat;
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

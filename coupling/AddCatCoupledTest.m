function tests = AddCatCoupledTest
    tests = functiontests(localfunctions);
end

function multipleTest(testCase)
    global BORROWING
    currBORROWING = BORROWING;
    for BORROWING = [0, 1]
        fprintf('Test BORROWING = %d\n', BORROWING);
        coupledTestRun(testCase);
        couplingTestRun(testCase);
        distributionTestRun(testCase);
    end
    BORROWING = currBORROWING;
end

function coupledTestRun(testCase)
    global BORROWING
    L = 6:10;
    n_i = length(L);
    n_j = 1e2;
    for i = 1:n_i
        [state_x, state_y] = unitTests.coupledStates(L(i), 1e-2);
        for j = 1:n_j
            [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
                    = AddCatCoupled(state_x, state_y);
            assertEqual(testCase, nstate_x.ncat, nstate_y.ncat);
            assertEqual(testCase, nstate_x.cat, nstate_y.cat);
            if BORROWING
                assertEqual(testCase, ...
                            [nstate_x.tree(find(nstate_x.cat)).catloc], ...
                            [nstate_y.tree(find(nstate_y.cat)).catloc]);
            end
            assertEqual(testCase, U_x, U_y);
            assertEqual(testCase, [OK_x, OK_y], [1, 1]);
            assertEqual(testCase, logq_x, logq_y);
            [nstate_z, U_z, OK_z, logq_z] = AddCat(state_x);
            if all(nstate_z.cat == nstate_x.cat)
                assertEqual(testCase, U_x, U_z);
                assertEqual(testCase, OK_x, OK_z);
                assertEqual(testCase, logq_x, logq_z);
            end
        end
    end
end

function couplingTestRun(testCase)
    L = 6:10;
    n_i = length(L);
    n_j = 1e4;
    [cObs, cExp] = deal(zeros(n_i, 1));
    for i = 1:n_i
        [state_x, state_y] = housekeptStates(L(i));
        for j = 1:n_j
            [nstate_x, nstate_y] = AddCatCoupled(state_x, state_y);
            if find(nstate_x.cat ~= state_x.cat) ...
                    == find(nstate_y.cat ~= state_y.cat)
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

function distributionTestRun(testCase)
    L = 6:10;
    n_i = length(L);
    n_j = 1e4;
    for i = 1:n_i
        [c_x, c_y, m_x, m_y] = deal(zeros(2 * L(i), 1));
        [state_x, state_y] = housekeptStates(L(i));
        for j = 1:n_j
            [c_nstate_x, c_nstate_y] = AddCatCoupled(state_x, state_y);

            ind_c_x = c_nstate_x.cat ~= state_x.cat;
            c_x(ind_c_x) = c_x(ind_c_x) + 1;

            ind_c_y = c_nstate_y.cat ~= state_y.cat;
            c_y(ind_c_y) = c_y(ind_c_y) + 1;

            m_nstate_x = AddCat(state_x);
            m_nstate_y = AddCat(state_y);

            ind_m_x = m_nstate_x.cat ~= state_x.cat;
            m_x(ind_m_x) = m_x(ind_m_x) + 1;

            ind_m_y = m_nstate_y.cat ~= state_y.cat;
            m_y(ind_m_y) = m_y(ind_m_y) + 1;
        end
        c_x = c_x / n_j;
        c_y = c_y / n_j;
        m_x = m_x / n_j;
        m_y = m_y / n_j;

        p_x = getBranchLengths(state_x) ./ state_x.length;
        p_y = getBranchLengths(state_y) ./ state_y.length;

        subplot(n_i, 2, 2 * i - 1);
        yyaxis left;
        plot([c_x, m_x, p_x'])
        yyaxis right;
        plot([c_x, m_x] - p_x');

        title(sprintf('L = %d : x', L(i)));
        if i == n_i
            xlabel('node')
            yyaxis left;
            ylabel('proportion selected');
        end

        subplot(n_i, 2, 2 * i);
        yyaxis left;
        plot([c_y, m_y, p_y'])
        yyaxis right;
        plot([c_y, m_y] - p_y');

        title(sprintf('L = %d : y', L(i)));
        if i == n_i
            legend('coupled', 'marginal', 'exact');
            yyaxis right;
            ylabel('difference from exact');
        end
    end
    v = input('Are the figures proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function [state_x, state_y] = housekeptStates(L)
    [state_x, state_y] = unitTests.housekeptStates(L, 1e-2);
    for i = 1:poissrnd(1)
        state_x = AddCat(state_x);
        state_y = AddCat(state_y);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT;
    MCMCCAT = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    clf;
end

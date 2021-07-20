function tests = DelCatCoupledTest
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
    M = length(L);
    N = 1e2;
    for m = 1:M
        [state_x, state_y] = unitTests.coupledStates(L(m), 1e-2);
        for c = 1:5
            [state_x, state_y] = AddCatCoupled(state_x, state_y);
        end
        for j = 1:N
            [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
                = DelCatCoupled(state_x, state_y);
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
            [nstate_z, U_z, OK_z, logq_z] = DelCat(state_x);
            if all(nstate_z.cat == nstate_x.cat)
                assertEqual(testCase, U_x, U_z);
                assertEqual(testCase, OK_x, OK_z);
                assertEqual(testCase, logq_x, logq_z);
            end
        end
    end
end

function couplingTestRun(testCase)
    global BORROWING;
    L = 5:10;
    M = length(L);
    N = 1e4;
    [cObs, cExp] = deal(zeros(M, 1));
    for m = 1:M
        [state_x, state_y] = housekeptStates(L(m));
        n_x = state_x.ncat;
        n_y = state_y.ncat;

        if BORROWING
            [locs_x, inds_x] = sampleCatIndex.getCats(state_x);
            [locs_y, inds_y] = sampleCatIndex.getCats(state_y);
        end

        for n = 1:N
            [nstate_x, nstate_y] = DelCatCoupled(state_x, state_y);
            if BORROWING
                [locn_x, ~] = sampleCatIndex.getCats(nstate_x);
                [locn_y, ~] = sampleCatIndex.getCats(nstate_y);

                ks_x = ~ismember(locs_x, locn_x);
                ks_y = ~ismember(locs_y, locn_y);

                i_x = inds_x(ks_x, 1);
                i_y = inds_y(ks_y, 1);

                % Same proposal if branch and location match
                if (i_x == i_y) && (locs_x(ks_x) == locs_y(ks_y))
                    cObs(m) = cObs(m) + 1;
                end
            else
                if (state_x.cat > nstate_x.cat) == (state_y.cat > nstate_y.cat)
                    cObs(m) = cObs(m) + 1;
                end
            end
        end
        if BORROWING
            for i = find(state_x.cat(:)' & state_y.cat(:)')
                n_xi = state_x.cat(i);
                n_yi = state_y.cat(i);
                l_x = state_x.tree(i).catloc;
                l_y = state_y.tree(i).catloc;
                cExp(m) = cExp(m) ...
                    + min(n_xi / n_x, n_yi / n_y) ...
                        * length(intersect(l_x, l_y)) / max(n_xi, n_yi);
            end
        else
            p_x = state_x.cat ./ n_x;
            p_y = state_y.cat ./ n_y;
            cExp(m) = sum(min(p_x, p_y));
        end
    end
    cObs = cObs / N;
    fprintf('Proportion of matching samples in each of %g trials\n', N);
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs - cExp));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function distributionTestRun(testCase)
    global BORROWING
    L = 6:10;
    M = length(L);
    N = 1e4;
    for m = 1:M
        [state_x, state_y] = housekeptStates(L(m));
        if BORROWING
            [c_x, m_x] = deal(zeros(state_x.ncat, 1));
            [c_y, m_y] = deal(zeros(state_y.ncat, 1));
            catloc_x = [state_x.tree.catloc];
            catloc_y = [state_y.tree.catloc];
        else
            [c_x, c_y, m_x, m_y] = deal(zeros(2 * L(m), 1));
        end
        for n = 1:N
            [c_nstate_x, c_nstate_y] = DelCatCoupled(state_x, state_y);
            m_nstate_x = DelCat(state_x);
            m_nstate_y = DelCat(state_y);

            if BORROWING
                [~, ind_c_x] = setdiff(catloc_x, [c_nstate_x.tree.catloc]);
                [~, ind_c_y] = setdiff(catloc_y, [c_nstate_y.tree.catloc]);
                [~, ind_m_x] = setdiff(catloc_x, [m_nstate_x.tree.catloc]);
                [~, ind_m_y] = setdiff(catloc_y, [m_nstate_y.tree.catloc]);
            else
                ind_c_x = c_nstate_x.cat ~= state_x.cat;
                ind_c_y = c_nstate_y.cat ~= state_y.cat;
                ind_m_x = m_nstate_x.cat ~= state_x.cat;
                ind_m_y = m_nstate_y.cat ~= state_y.cat;
            end
            c_x(ind_c_x) = c_x(ind_c_x) + 1;
            c_y(ind_c_y) = c_y(ind_c_y) + 1;
            m_x(ind_m_x) = m_x(ind_m_x) + 1;
            m_y(ind_m_y) = m_y(ind_m_y) + 1;
        end
        c_x = c_x / N;
        c_y = c_y / N;
        m_x = m_x / N;
        m_y = m_y / N;

        if BORROWING
            p_x = ones(state_x.ncat, 1) / state_x.ncat;
            p_y = ones(state_y.ncat, 1) / state_y.ncat;
        else
            p_x = state_x.cat ./ state_x.ncat;
            p_y = state_y.cat ./ state_y.ncat;
        end

        i_x = find(p_x > 0);
        i_y = find(p_y > 0);

        subplot(M, 2, 2 * m - 1);
        plot(i_x, [c_x(i_x), m_x(i_x)] - p_x(i_x), 'x');
        xticks(unique(round(get(gca, 'xTick'))));

        title(sprintf('L = %d : x : %d cats', L(m), state_x.ncat));
        if m == M
            if BORROWING
                xlabel('catastrophe');
            else
                xlabel('catastrophe branch')
            end
            ylabel('difference from exact');
        end

        subplot(M, 2, 2 * m);
        plot(i_y, [c_y(i_y), m_y(i_y)] - p_y(i_y), 'x');
        xticks(unique(round(get(gca, 'xTick'))));

        title(sprintf('L = %d : y : %d cats', L(m), state_y.ncat));
        if m == M
            if BORROWING
                xlabel('catastrophe');
            else
                xlabel('catastrophe branch')
            end
            legend('coupled', 'marginal');
        end
    end
    v = input('Are the figures proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function [state_x, state_y] = housekeptStates(L)
    [state_x, state_y] = unitTests.housekeptStates(L, 1e-2);
    for c = 1:(3 + poissrnd(1))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    for c = 1:(3 + poissrnd(1))
        state_x = AddCat(state_x);
    end
    for c = 1:(3 + poissrnd(1))
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
    close;
end

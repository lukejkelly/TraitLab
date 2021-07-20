function tests = MoveCatCoupledTest
    tests = functiontests(localfunctions);
end

function multipleTest(testCase)
    global BORROWING
    currBORROWING = BORROWING;
    for BORROWING = [0, 1]
        fprintf('Test BORROWING = %d\n', BORROWING);
        coupledTestRun(testCase);
        couplingTestRun(testCase);
    end
    BORROWING = currBORROWING;
end

function coupledTestRun(testCase)
    global BORROWING
    L = 4:10;
    M = length(L);
    N = 1e2;
    for m = 1:M
        [state_x, state_y] = unitTests.coupledStates(L(m), 1e-2);
        for c = 1:5
            [state_x, state_y] = AddCatCoupled(state_x, state_y);
        end
        for n = 1:N
            [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
                    = MoveCatCoupled(state_x, state_y);
            assertEqual(testCase, nstate_x.ncat, nstate_y.ncat);
            assertEqual(testCase, nstate_x.cat, nstate_y.cat);
            if BORROWING
                cellfun(@(x, y) assertEqual(testCase, x, y), ...
                        {nstate_x.tree.catloc}, {nstate_y.tree.catloc});
            end
            assertEqual(testCase, U_x, U_y);
            assertEqual(testCase, [OK_x, OK_y], [1, 1]);
            assertEqual(testCase, logq_x, logq_y);
            [nstate_z, U_z, OK_z, logq_z] = MoveCat(state_x);
            if all(nstate_z.cat == nstate_x.cat)
                assertEqual(testCase, U_x, U_z);
                assertEqual(testCase, OK_x, OK_z);
                assertEqual(testCase, logq_x, logq_z);
            end
        end
    end
end

function couplingTestRun(testCase)
    global BORROWING ROOT;
    L = 4:10;
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
            [nstate_x, nstate_y] = MoveCatCoupled(state_x, state_y);
            if BORROWING
                [locn_x, indn_x] = sampleCatIndex.getCats(nstate_x);
                [locn_y, indn_y] = sampleCatIndex.getCats(nstate_y);

                ks_x = ~ismember(locs_x, locn_x);
                ks_y = ~ismember(locs_y, locn_y);

                i_x = inds_x(ks_x, 1);
                i_y = inds_y(ks_y, 1);

                kn_x = ~ismember(locn_x, locs_x);
                kn_y = ~ismember(locn_y, locs_y);

                j_x = indn_x(kn_x, 1);
                j_y = indn_y(kn_y, 1);

                % If branch (source and destination) and location match
                if (i_x == i_y) && (locs_x(ks_x) == locs_y(ks_y)) && (j_x == j_y)
                    cObs(m) = cObs(m) + 1;
                end
            else
                if all((nstate_x.cat - state_x.cat) ...
                        == (nstate_y.cat - state_y.cat))
                   cObs(m) = cObs(m) + 1;
                end
            end
        end
        root = state_x.root;
        if BORROWING
            for i = find(state_x.cat(:) & state_y.cat(:))'
                n_xi = state_x.cat(i);
                n_yi = state_y.cat(i);
                o_b = min(n_xi / n_x, n_yi / n_y);

                l_xi = locs_x(inds_x(:, 1) == i);
                l_yi = locs_y(inds_y(:, 1) == i);
                o_l = length(intersect(l_xi, l_yi)) / max(n_xi, n_yi);

                [u_x, v_x] = GetLegalCoupled.getPoss(state_x.tree, i, root);
                [u_y, v_y] = GetLegalCoupled.getPoss(state_y.tree, i, root);
                n_u = length(intersect(u_x, u_y));
                n_v = max(v_x, v_y);
                o_d = n_u / n_v;

                cExp(m) = cExp(m) + o_b * o_l * o_d;
            end
        else
            if state_x.ncat == 0 || state_y.ncat == 0
                cExp(m) = 1;
            else
                p_x = state_x.cat ./ state_x.ncat;
                p_y = state_y.cat ./ state_y.ncat;
                for k = 1:(2 * L(m))
                    if state_x.tree(k).type < ROOT
                        [u_x, v_x] = GetLegalCoupled.getPoss(state_x.tree, k, root);
                        [u_y, v_y] = GetLegalCoupled.getPoss(state_y.tree, k, root);
                        n_u = length(intersect(u_x, u_y));
                        n_v = max(v_x, v_y);
                        cExp(m) = cExp(m) + min(p_x(k), p_y(k)) * n_u / n_v;
                    end
                end
            end
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
end

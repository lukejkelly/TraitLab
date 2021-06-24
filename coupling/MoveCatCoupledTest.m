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
    n_i = length(L);
    n_j = 1e2;
    for i = 1:n_i
        [state_x, state_y] = unitTests.coupledStates(L(i), 1e-2);
        for c = 1:5
            [state_x, state_y] = AddCatCoupled(state_x, state_y);
        end
        for j = 1:n_j
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
    n_i = length(L);
    n_j = 1e4;
    [cObs, cExp] = deal(zeros(n_i, 1));
    for i = 1:n_i
        [state_x, state_y] = housekeptStates(L(i));
        if BORROWING
            catloc_x = {state_x.tree.catloc};
            catloc_y = {state_y.tree.catloc};
        end
        for j = 1:n_j
            [nstate_x, nstate_y] = MoveCatCoupled(state_x, state_y);
            if BORROWING
                cat_x = cellfun(@setdiff, catloc_x, {nstate_x.tree.catloc}, ...
                                'UniformOutput', false);
                cat_y = cellfun(@setdiff, catloc_y, {nstate_y.tree.catloc}, ...
                                'UniformOutput', false);
                cat = ismembertol([cat_x{:}], [cat_y{:}]);
                dest_x = find(nstate_x.cat > state_x.cat);
                dest_y = find(nstate_y.cat > state_y.cat);
                dest = dest_x == dest_y;
                if cat && dest
                    cObs(i) = cObs(i) + 1;
                end
            else
                if all((nstate_x.cat - state_x.cat) ...
                       == (nstate_y.cat - state_y.cat))
                   cObs(i) = cObs(i) + 1;
               end
            end
        end
        root = state_x.root;
        if BORROWING
            i_c = intersect([catloc_x{:}], [catloc_y{:}]);
            l_i = length(i_c);
            n_c = max(state_x.ncat, state_y.ncat);
            if l_i > 0
                s_x = state_x.tree;
                s_y = state_y.tree;
                for k = 1:l_i
                    k_x = find(cellfun(@(x) ismembertol(i_c(k), x), {s_x.catloc}));
                    k_y = find(cellfun(@(x) ismembertol(i_c(k), x), {s_y.catloc}));
                    [u_x, v_x] = GetLegalCoupled.getPoss(s_x, k_x, root);
                    [u_y, v_y] = GetLegalCoupled.getPoss(s_y, k_y, root);
                    n_u = length(intersect(u_x, u_y));
                    n_v = max(v_x, v_y);
                    cExp(i) = cExp(i) + (1 / n_c) * n_u / n_v;
                end
            end
        else
            if state_x.ncat == 0 || state_y.ncat == 0
                cExp(i) = 1;
            else
                p_x = state_x.cat ./ state_x.ncat;
                p_y = state_y.cat ./ state_y.ncat;
                for k = 1:(2 * L(i))
                    if state_x.tree(k).type < ROOT
                        [u_x, v_x] = GetLegalCoupled.getPoss(state_x.tree, k, root);
                        [u_y, v_y] = GetLegalCoupled.getPoss(state_y.tree, k, root);
                        n_u = length(intersect(u_x, u_y));
                        n_v = max(v_x, v_y);
                        cExp(i) = cExp(i) + min(p_x(k), p_y(k)) * n_u / n_v;
                    end
                end
            end
        end
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
    for c = 1:(1 + poissrnd(1))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    for c = 1:(1 + poissrnd(1))
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
end

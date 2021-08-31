function tests = BcatsCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global BORROWING NARROW WIDE
    M = 10;
    N = 10;
    for BORROWING = [0, 1]
        for mt = [NARROW, WIDE]
            for addClades = [0, 1]
                for m = 1:M
                    L = 8 + floor(rand * 3);
                    [state_x, state_y, theta, prior] = coupledStates(L, addClades);
                    for n = 1:N
                        newage = [];
                        while isempty(newage)
                            [i, j, k, newage, ~, ~, ~, ~] = Bchoose(...
                                state_x, mt, theta, prior);
                        end
                        [ncat_x, ncat_y, cat_x, cat_y, loc_x, loc_y, ...
                         logq_x, logq_y] = BcatsCoupled(...
                            state_x, state_y, i, j, j, k, k, newage, newage);

                        assertEqual(testCase, ncat_x, ncat_y);
                        assertEqual(testCase, cat_x, cat_y);
                        assertEqual(testCase, loc_x, loc_y);
                        assertEqual(testCase, logq_x, logq_y);
                    end
                end
            end
        end
    end
end

function couplingTest(testCase)
    % For a valid SPR move, check if distributions overlap correctly
    global BORROWING WIDE
    M = 8;
    N = 1e4;
    L = 10;
    [gObs, gExp] = deal(nan(M, 1));
    BORROWING = 0;
    mt = WIDE;
    addClades = 0;
    for m = 1:M
        [state_x, state_y, theta, prior] = housekeptStates(L, addClades);
        r = state_x.root;
        [newage_x, newage_y] = deal([]);
        while isempty(newage_x) || isempty(newage_y)
            [i, j_x, j_y, k_x, k_y, newage_x, newage_y, ~, ~, ~, ~, ~, ~, ~, ~] ...
                = BchooseCoupled(state_x, state_y, mt, theta, prior);
        end
        fObs = 0;
        mObs = 0;
        for n = 1:N
            [~, ~, cat_x, cat_y, ~, ~, ~, ~] = BcatsCoupled(...
                state_x, state_y, i, j_x, j_y, k_x, k_y, newage_x, newage_y);
            if cat_x.j == cat_y.j
                fObs = fObs + 1;
            end
            mj = max(cat_x.j, cat_y.j);
            if mj > mObs
                mObs = mj;
            end
        end
        gObs(m) = fObs / N;

        s_x = state_x.tree;
        dj_x = newage_x - s_x(j_x).time;
        wj_x = dj_x / (s_x(k_x).time - s_x(j_x).time);
        if j_x == r
            f_x = @(z) exp(Bcats.logPriorCounts(z, dj_x, state_x.rho));
        else
            f_x = @(z) binopdf(z, state_x.cat(j_x), wj_x);
        end

        s_y = state_y.tree;
        dj_y = newage_y - s_y(j_y).time;
        wj_y = dj_y / (s_y(k_y).time - s_y(j_y).time);
        if j_y == r
            f_y = @(z) exp(Bcats.logPriorCounts(z, dj_y, state_y.rho));
        else
            f_y = @(z) binopdf(z, state_y.cat(j_y), wj_y);
        end

        gExp(m) = sum(arrayfun(@(z) min(f_x(z), f_y(z)), 0:mObs));
    end

    fprintf('Proportion of matching samples in each of %g trials\n', N);
    fprintf('Repeated for %g 10-leaf trees\n\n', M);

    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', gObs(:, 1)));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', gExp(:, 1)));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', gObs(:, 1) - gExp(:, 1)));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function [state_x, state_y, theta, prior] = coupledStates(L, addClades)
    global ROOT
    theta = 1e-3 * rand;
    [state_x, state_y] = unitTests.coupledStates(L, theta);
    if addClades
        s = state_x.tree;

        clade = synthclades(s, 1, 2, 1 - rand^3);
        rootmax = 2 * s([s.type] == ROOT).time;
        prior = unitTests.clade2prior(clade, rootmax);

        state_x.tree = treeclades(state_x.tree, prior.clade);
        state_x = UpdateClades(state_x, [state_x.leaves, state_x.nodes], ...
                               size(prior.clade, 2));

        state_y.tree = treeclades(state_y.tree, prior.clade);
        state_y = UpdateClades(state_y, [state_y.leaves, state_y.nodes], ...
                               size(prior.clade, 2));
    else
        prior = struct('isclade', 0);
    end
    for c = 1:(3 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    state_x.rho = state_x.ncat / state_x.length;
    state_y.rho = state_y.ncat / state_y.length;
end

function [state_x, state_y, theta, prior] = housekeptStates(L, addClades)
    global ROOT
    theta = rand * 1e-3;
    if addClades
        state = unitTests.dummyState(ExpTree(L, theta));
        s = state.tree;
        clade = synthclades(s, 1, 2, 1 - rand^3);
        rootmax = (1 + rand) * s([s.type] == ROOT).time;
        prior = unitTests.clade2prior(clade, rootmax);

        s_x = RandCladeTree(clade, {s(state.leaves).Name}, rootmax, 0);
        s_x = treeclades(s_x, prior.clade);
        state_x = unitTests.dummyState(s_x);
        state_x = UpdateClades(state_x, [state_x.leaves, state_x.nodes], ...
                               size(prior.clade, 2));

        s_y = RandCladeTree(clade, {s(state.leaves).Name}, rootmax, 0);
        s_y = treeclades(s_y, prior.clade);
        state_y = unitTests.dummyState(s_y);
        state_y = UpdateClades(state_y, [state_y.leaves, state_y.nodes], ...
                               size(prior.clade, 2));
        state_y = housekeeping(state_x, state_y);
    else
       [state_x, state_y] = unitTests.housekeptStates(L, 5e-3 * rand);
       prior = struct('isclade', 0);
    end
    for c = 1:(3 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    for c = 1:(3 + poissrnd(3))
        state_x = AddCat(state_x);
    end
    for c = 1:(3 + poissrnd(3))
        state_y = AddCat(state_y);
    end
    state_x.rho = state_x.ncat / state_x.length;
    state_y.rho = state_y.ncat / state_y.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    global MCMCCAT VARYRHO
    unitTests.setupOnce(testCase);
    testCase.TestData.VARYRHO = VARYRHO;
    MCMCCAT = 1;
    VARYRHO = 1;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

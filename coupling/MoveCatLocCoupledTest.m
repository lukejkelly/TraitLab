function tests = MoveCatLocCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    M = 1e3;
    [ncat, cat, i, U, OK, logq] = deal(nan(M, 1));
    for m = 1:M
        L = 8 + ceil(rand * 8);
        [state_x, state_y] = coupledStates(L);
        [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
            = MoveCatLocCoupled(state_x, state_y);

        ncat(m) = nstate_x.ncat == nstate_y.ncat;
        cat(m) = all(nstate_x.cat == nstate_y.cat);

        [locs_x, inds_x] = sampleCatIndex.getCats(state_x);
        [locs_y, inds_y] = sampleCatIndex.getCats(state_y);

        [locn_x, indn_x] = sampleCatIndex.getCats(nstate_x);
        [locn_y, indn_y] = sampleCatIndex.getCats(nstate_y);

        i(m) = all(locs_x == locs_y) && all(inds_x == inds_y, 'all') ...
               && all(locn_x == locn_y) && all(indn_x == indn_y, 'all');

        U(m) = all(U_x == U_y);
        OK(m) = OK_x == 1 && OK_x == OK_y;
        logq(m) = logq_x == 0 && logq_x == logq_y;
    end
    assertTrue(testCase, all(ncat));
    assertTrue(testCase, all(cat));
    assertTrue(testCase, all(i));
    assertTrue(testCase, all(U));
    assertTrue(testCase, all(OK));
    assertTrue(testCase, all(logq));
end

function couplingTest(testCase)
    % Check marginal and coupling selection distributions
    M = 8;
    N = 1e4;
    [i_x, i_y, catloc_new] = deal(nan(M, N));
    fprintf('Observed and expected coupling proportions\n');
    fprintf('%-4s%-8s%-8s%-8s\n', 'L', 'Obs', 'Exp', 'Diff');
    for m = 1:M
        L = 8 + ceil(rand * 8);
        [state_x, state_y] = housekeptStates(L);
        n_x = state_x.ncat;
        n_y = state_y.ncat;
        [locs_x, inds_x] = sampleCatIndex.getCats(state_x);
        [locs_y, inds_y] = sampleCatIndex.getCats(state_y);

        fObs_x = zeros(1, state_x.ncat);
        fObs_y = zeros(1, state_y.ncat);
        fObs_c = 0;

        for n = 1:N
            [nstate_x, nstate_y, ~, ~, ~, ~, ~, ~] ...
                = MoveCatLocCoupled(state_x, state_y);

            [locn_x, indn_x] = sampleCatIndex.getCats(nstate_x);
            [locn_y, indn_y] = sampleCatIndex.getCats(nstate_y);

            ks_x = ~ismember(locs_x, locn_x);
            ks_y = ~ismember(locs_y, locn_y);

            fObs_x(ks_x) = fObs_x(ks_x) + 1;
            fObs_y(ks_y) = fObs_y(ks_y) + 1;

            if inds_x(ks_x, 1) == inds_y(ks_y) ...
                    && ismembertol(locs_x(ks_x), locs_y(ks_y))
                fObs_c = fObs_c + 1;
            end

            kn_x = ~ismember(locn_x, locs_x);
            kn_y = ~ismember(locn_y, locs_y);

            i_x(m, n) = indn_x(kn_x, 1) == inds_x(ks_x, 1);
            i_y(m, n) = indn_y(kn_y, 1) == inds_y(ks_y, 1);
            catloc_new(m, n) = locn_x(kn_x) == locn_y(kn_y);
        end

        gObs_x = cumsum(fObs_x) / N;
        gObs_y = cumsum(fObs_y) / N;
        gObs_c = fObs_c / N;

        gExp_x = (1:state_x.ncat) / state_x.ncat;
        gExp_y = (1:state_y.ncat) / state_y.ncat;
        gExp_c = 0;
        for i = find(state_x.cat & state_y.cat)'
            l_xi = locs_x(inds_x(:, 1) == i);
            l_yi = locs_y(inds_y(:, 1) == i);
            n_xi = state_x.cat(i);
            n_yi = state_y.cat(i);

            gExp_c = gExp_c + length(intersect(l_xi, l_yi)) / max(n_xi, n_yi) ...
                              * min(n_xi / n_x, n_yi / n_y);
        end

        nexttile;
        plot(1:state_x.ncat, gObs_x, '*', 1:state_x.ncat, gExp_x, 'o')
        title(sprintf('L = %d : x', L));

        nexttile;
        plot(1:state_y.ncat, gObs_y, '*', 1:state_y.ncat, gExp_y, 'o')
        title(sprintf('L = %d : y', L));

        fprintf('%-4d%-8.4g%-8.4g%-8.4g\n', L, gObs_c, gExp_c, gObs_c - gExp_c);
    end
    assertTrue(testCase, all(catloc_new, 'all'));
    fprintf('\nFigure displays marginal distributions\n');
    v = input('Do these proportions match? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function [state_x, state_y] = coupledStates(L)
    [state_x, state_y] = unitTests.coupledStates(L, 5e-3 * rand);
    for c = 1:(1 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
end

function [state_x, state_y] = housekeptStates(L)
    [state_x, state_y] = unitTests.housekeptStates(L, 5e-3 * rand);
    for c = 1:(1 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    for c = 1:(1 + poissrnd(3))
        state_x = AddCat(state_x);
    end
    for c = 1:(1 + poissrnd(3))
        state_y = AddCat(state_y);
    end
    for c = 1:(1 + poissrnd(3))
        [state_x, state_y] = MoveCatCoupled(state_x, state_y);
    end
    for c = 1:(1 + poissrnd(3))
        state_x = MoveCat(state_x);
    end
    for c = 1:(1 + poissrnd(3))
        state_y = MoveCat(state_y);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 1;
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

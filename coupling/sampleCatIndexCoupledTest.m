function tests = sampleCatIndexCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    M = 1e2;
    N = 1e2;
    for m = 1:M
        L = 5 + poissrnd(2);
        [state_x, state_y] = coupledStates(L);
        c = find(state_x.cat);
        C = length(c);
        [j_x, j_y] = deal(nan(C, N));
        for h = 1:C
            i = c(h);
            for n = 1:N
                [j_x(h, n), j_y(h, n)] ...
                    = sampleCatIndexCoupled(state_x, state_y, i, i);
            end
        end
        assertEqual(testCase, j_y, j_y);
    end
end

function housekeptTest(testCase)
    M = 4;
    N = 1e4;
    [cObs, cExp] = deal(zeros(M, 1));
    for m = 1:M
        L = 5 + poissrnd(2);
        [state_x, state_y] = housekeptStates(L);
        n_x = state_x.ncat;
        n_y = state_y.ncat;

        [locs_x, inds_x] = sampleCatIndex.getCats(state_x);
        [locs_y, inds_y] = sampleCatIndex.getCats(state_y);

        fObs_x = zeros(n_x, 1);
        fObs_y = zeros(n_y, 1);

        for n = 1:N
            [i_x, i_y] = sampleBranchByCatCountCoupled(state_x, state_y);
            [j_x, j_y] = sampleCatIndexCoupled(state_x, state_y, i_x, i_y);

            k_x = i_x == inds_x(:, 1) & j_x == inds_x(:, 2);
            k_y = i_y == inds_y(:, 1) & j_y == inds_y(:, 2);

            fObs_x(k_x) = fObs_x(k_x) + 1;
            fObs_y(k_y) = fObs_y(k_y) + 1;

            if i_x == i_y && locs_x(k_x) == locs_y(k_y)
                cObs(m) = cObs(m) + 1;
            end
        end
        gObs_x = cumsum(fObs_x / N);
        gObs_y = cumsum(fObs_y / N);

        gExp_x = (1:n_x) / n_x;
        gExp_y = (1:n_y) / n_y;

        for i = find(state_x.cat & state_y.cat)'
            l_xi = locs_x(inds_x(:, 1) == i);
            l_yi = locs_y(inds_y(:, 1) == i);
            n_xi = state_x.cat(i);
            n_yi = state_y.cat(i);

            cExp(m) = cExp(m) + length(intersect(l_xi, l_yi)) / max(n_xi, n_yi) ...
                                * min(n_xi / n_x, n_yi / n_y);
        end

        subplot(M, 2, 2 * m - 1);
        plot(1:n_x, gObs_x, 'o', 1:n_x, gExp_x, 'x');
        xlabel('Catastrophe');
        ylabel('CDF');
        title(sprintf('x: %d leaves, %d cats', L, n_x));

        subplot(M, 2, 2 * m);
        plot(1:n_y, gObs_y, 'o', 1:n_y, gExp_y, 'x');
        xlabel('Catastrophe');
        title(sprintf('y: %d leaves, %d cats', L, n_y));
    end
    legend('ECDF', 'CDF', 'Location', 'southeast');
    fprintf('Observed and expected proportions from %g samples\n', N);
    v1 = input('Are these marginal CDFs okay? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);

    cObs = cObs / N;
    fprintf('Proportion of matching samples in each of %g trials\n', N);
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs - cExp));
    v2 = input('Are these coupling proportions okay? Reply 1 for yes... ');
    assertTrue(testCase, v2 == 1);
end

function [state_x, state_y] = housekeptStates(L)
    [state_x, state_y] = unitTests.housekeptStates(L, 1e-2 * rand);
    for c = 1:(3 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    for c = 1:(3 + poissrnd(3))
        state_x = AddCat(state_x);
    end
    for c = 1:(3 + poissrnd(3))
        state_y = AddCat(state_y);
    end
end

function [state_x, state_y] = coupledStates(L)
    [state_x, state_y] = unitTests.coupledStates(L, 1e-2 * rand);
    for c = 1:(3 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
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

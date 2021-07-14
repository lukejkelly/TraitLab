function tests = sampleCatIndexCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    n_i = 1e2;
    n_j = 1e2;
    for i = 1:n_i
        L = 5 + poissrnd(2);
        [state_x, state_y] = coupledStates(L);
        [i_x, i_y, j_x, j_y] = deal(nan(n_j, 1));
        for j = 1:n_j
            [i_x(j), i_y(j), j_x(j), j_y(j)] ...
                = sampleCatIndexCoupled(state_x, state_y);
        end
        assertEqual(testCase, i_x, i_y);
        assertEqual(testCase, j_y, j_y);
    end
end


function housekeptTest(testCase)
    n_i = 4;
    n_j = 1e5;
    [cObs, cExp] = deal(zeros(n_i, 1));
    for i = 1:n_i
        L = 5 + poissrnd(2);
        [state_x, state_y] = housekeptStates(L);

        [locs_x, inds_x] = getCats(state_x);
        [locs_y, inds_y] = getCats(state_y);

        f_x = zeros(state_x.ncat, 1);
        f_y = zeros(state_y.ncat, 1);

        for j = 1:n_j
            [i_x, i_y, j_x, j_y] = sampleCatIndexCoupled(state_x, state_y);
            k_x = find(i_x == inds_x(:, 1) & j_x == inds_x(:, 2));
            k_y = find(i_y == inds_y(:, 1) & j_y == inds_y(:, 2));

            f_x(k_x) = f_x(k_x) + 1;
            f_y(k_y) = f_y(k_y) + 1;

            if locs_x(k_x) == locs_y(k_y)
                cObs(i) = cObs(i) + 1;
            end
        end

        cExp(i) = length(intersect(locs_x, locs_y)) ...
                  / max(state_x.ncat, state_y.ncat);

        subplot(n_i, 2, 2 * i - 1);
        plot(1:state_x.ncat, cumsum(f_x) / n_j, 'o', ...
             1:state_x.ncat, cumsum(ones(1, state_x.ncat)) / state_x.ncat, 'x');
        legend('ECDF', 'CDF');
        xlabel('Catastrophe');
        ylabel('CDF');
        title(sprintf('x: %d leaves, %d cats', L, state_x.ncat));

        subplot(n_i, 2, 2 * i);
        plot(1:state_y.ncat, cumsum(f_y) / n_j, 'o', ...
             1:state_y.ncat, cumsum(ones(1, state_y.ncat)) / state_y.ncat, 'x');
        legend('ECDF', 'CDF');
        xlabel('Catastrophe');
        ylabel('CDF');
        title(sprintf('y: %d leaves, %d cats', L, state_y.ncat));
    end
    fprintf('Observed and expected proportions from %g samples\n', n_j);
    v1 = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);

    cObs = cObs / n_j;
    fprintf('Proportion of matching samples in each of %g trials\n', n_j);
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs - cExp));
    v2 = input('Are these proportions the same? Reply 1 for yes... ');
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

function [locs, inds] = getCats(state)
    locs = nan(state.ncat, 1);
    inds = nan(state.ncat, 2);
    k = 1;
    for j = find(state.cat(:)')
        jInds = k:(k + state.cat(j) - 1);
        locs(jInds) = state.tree(j).catloc;
        inds(jInds, 1) = j;
        inds(jInds, 2) = 1:state.cat(j);
        k = k + state.cat(j);
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

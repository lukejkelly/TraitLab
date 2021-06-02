function tests = RscaleSubTreeCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    % Only checking variation terms
    rangeL = 6:10;
    theta = 1e-2;
    n_i = length(rangeL);
    n_j = 1e4;
    del = 0.5;
    deldel = 1.5;
    a = del;
    b = del + deldel;
    [o1, o2, c1, c2] = deal(nan(n_i, n_j));
    for i = 1:n_i
        L = rangeL(i);
        % Housekept
        [state_x1, state_y1] = RscaleSubTreeCoupled.unitTests.housekeptStates(L, theta);
        % Coupled
        [state_x2, state_y2] = RscaleSubTreeCoupled.unitTests.coupledStates(L, theta);
        for j = 1:n_j
            % Housekept
            [nstate_x1, nstate_y1, ~, ~, ~, ~, ~, ~] = RscaleSubTreeCoupled(...
                state_x1, state_y1, del, deldel);
            i_x1 = findSubTreeRoot(state_x1, nstate_x1);
            i_y1 = findSubTreeRoot(state_y1, nstate_y1);

            if rand < 1e-2
                rho_x1 = getRho(i_x1, state_x1, nstate_x1);
                rho_y1 = getRho(i_y1, state_y1, nstate_y1);
                checkRange(testCase, rho_x1, rho_y1, a, b);
                assertNotEqual(testCase, rho_x1, rho_y1);
            end

            c1(i, j) = ismembertol(nstate_x1.tree(i_x1).time, nstate_y1.tree(i_y1).time);
            o1(i, j) = getOverlap(state_x1, state_y1, a, b, i_x1, i_y1);

            % Already coupled
            [nstate_x2, nstate_y2, ~, ~, ~, ~, ~, ~] = RscaleSubTreeCoupled(...
                state_x2, state_y2, del, deldel);
            i2 = findSubTreeRoot(state_x2, nstate_x2);
            assertEqual(testCase, i2, findSubTreeRoot(state_y2, nstate_y2));

            if rand < 1e-2
                rho_x2 = getRho(i2, state_x2, nstate_x2);
                rho_y2 = getRho(i2, state_y2, nstate_y2);
                checkRange(testCase, rho_x2, rho_y2, a, b);
                assertEqual(testCase, rho_x2, rho_y2);
            end

            c2(i, j) = ismembertol(nstate_x2.tree(i2).time, nstate_y2.tree(i2).time);
            o2(i, j) = getOverlap(state_x2, state_y2, a, b, i2, i2);
        end
    end

    % Proportion of matching samples
    fprintf('Proportion of matching samples after %g trials\n', n_i * n_j);
    fprintf('Housekept trees\n');
    fprintf('  Theoretical          = %g\n', mean(o1(:)));
    fprintf('  RscaleSubTreeCoupled = %g\n', mean(c1(:)));
    fprintf('Coupled trees\n');
    fprintf('  Theoretical          = %g\n', mean(o2(:)));
    fprintf('  RscaleSubTreeCoupled = %g\n', mean(c2(:)));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function marginalTest(testCase)
    rangeL = 7:11;
    theta = 1e-2;
    n_i = length(rangeL);
    n_j = 1e4;
    del = 0.5;
    deldel = 1.5;
    a = del;
    b = del + deldel;
    [xc, yc, xm, ym] = deal(nan(n_i, n_j));
    for i = 1:n_i
        L = rangeL(i);
        [state_x, state_y] = RscaleSubTreeCoupled.unitTests.housekeptStates(L, theta);
        for j = 1:n_j
            % Coupled
            [nstate_xc, nstate_yc, ~, ~, ~, ~, ~, ~] = RscaleSubTreeCoupled(...
                state_x, state_y, del, deldel);

            i_xc = findSubTreeRoot(state_x, nstate_xc);
            i_yc = findSubTreeRoot(state_y, nstate_yc);
            rho_xc = getRho(i_xc, state_x, nstate_xc);
            rho_yc = getRho(i_yc, state_y, nstate_yc);

            if rand < 1e-2
                checkRange(testCase, rho_xc, rho_yc, a, b);
            end

            xc(i, j) = rho_xc;
            yc(i, j) = rho_yc;

            % Marginal
            [nstate_xm, ~, ~, ~, ~] = RscaleSubTree(state_x, del, deldel);
            [nstate_ym, ~, ~, ~, ~] = RscaleSubTree(state_y, del, deldel);

            i_xm = findSubTreeRoot(state_x, nstate_xm);
            i_ym = findSubTreeRoot(state_y, nstate_ym);
            rho_xm = getRho(i_xm, state_x, nstate_xm);
            rho_ym = getRho(i_ym, state_y, nstate_ym);

            if rand < 1e-2
                checkRange(testCase, rho_xm, rho_ym, a, b);
            end

            xm(i, j) = rho_xm;
            ym(i, j) = rho_ym;
        end
    end

    % Empirical CDFs
    [ecdf_xc, x_xc] = ecdf(xc(:));
    [ecdf_yc, x_yc] = ecdf(yc(:));

    [ecdf_xm, x_xm] = ecdf(xm(:));
    [ecdf_ym, x_ym] = ecdf(ym(:));

    % Actual CDFs
    f_cdf = @(t) (min(max(t, a), b) - a) / (b - a);

    cdf_xc = arrayfun(f_cdf, x_xc(:));
    cdf_yc = arrayfun(f_cdf, x_yc(:));
    cdf_xm = arrayfun(f_cdf, x_xm(:));
    cdf_ym = arrayfun(f_cdf, x_ym(:));

    % Make figures
    subplot(2, 1, 1);
    plot([x_xc, x_xm], [ecdf_xc - cdf_xc, ecdf_xm - cdf_xm]);
    title('x');

    subplot(2, 1, 2);
    plot([x_yc, x_ym], [ecdf_yc - cdf_yc, ecdf_ym - cdf_ym]);
    title('y');

    for k = 1:2
        subplot(2, 1, k);
        xlabel('New root age');
        ylabel('Difference from exact CDF');
        legend('coupled', 'marginal');
    end

    fprintf('Difference from CDF for coupled and marginal samples\n');
    fprintf('%g draws for each of %g pairs of trees\n', n_j, n_i);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

% Helper functions
function ol = getOverlap(state_x, state_y, a, b, i_x, i_y)
    global LEAF;
    nu_x = progeny(state_x.tree, i_x, LEAF);
    nu_y = progeny(state_y.tree, i_y, LEAF);

    x_c = state_x.tree(i_x).time;
    y_c = state_y.tree(i_y).time;

    t0_x = min([state_x.tree(nu_x(1, :)).time]);
    t0_y = min([state_y.tree(nu_y(1, :)).time]);

    a_p = t0_x + a * (x_c - t0_x);
    b_p = t0_x + b * (x_c - t0_x);

    a_q = t0_y + a * (y_c - t0_y);
    b_q = t0_y + b * (y_c - t0_y);

    max_a = max(a_p, a_q);
    min_b = min(b_p, b_q);

    if min_b < max_a
        ol = 0;
    else
        ol = (min_b - max_a) / ((b - a) * max(x_c - t0_x, y_c - t0_y));
    end
end

function rho = getRho(i, state, nstate)
    global LEAF;
    nu = progeny(state.tree, i, LEAF);
    t0 = min([state.tree(nu(1, :)).time]);
    rho = (nstate.tree(i).time - t0) / (state.tree(i).time - t0);
end

function i = findSubTreeRoot(state, nstate)
    nInds = state.nodes;
    sInds = nInds([state.tree(nInds).time] ~= [nstate.tree(nInds).time]);
    [~, j] = max([state.tree(sInds).time]);
    i = sInds(j);
end

function checkRange(testCase, var_x, var_y, a, b)
    assertGreaterThanOrEqual(testCase, [var_x, var_y], a);
    assertLessThanOrEqual(testCase, [var_x, var_y], b);
end


% Setup and teardown functions
function setupOnce(testCase)
    RscaleSubTreeCoupled.unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    RscaleSubTreeCoupled.unitTests.teardownOnce(testCase);
end

function tests = RscaleTopTreeCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    % Only checking variation terms
    load('+RscaleTopTreeCoupledTest/cladePrior.mat', 'prior');
    rangeL = repelem(6:10, 2);
    n_i = length(rangeL);
    n_j = 2e3;
    del = 0.5;
    deldel = 1.5;
    a = del;
    b = del + deldel;
    [o1, o2, c1, c2] = deal(nan(n_i, n_j));
    for i = 1:n_i
        L = rangeL(i);
        % Housekept
        [state_x1, state_y1] = RscaleTopTreeCoupledTest.housekeptStates(L);
        r1 = state_x1.root;
        % Coupled
        [state_x2, state_y2] = RscaleTopTreeCoupledTest.coupledStates(L);
        r2 = state_x2.root;
        for j = 1:n_j
            % Housekept
            [nstate_x1, nstate_y1, ~, ~, ~, ~, ~, ~] = RscaleTopTreeCoupled(...
                state_x1, state_y1, prior, del, deldel);

            rho_x1 = getRho(state_x1, nstate_x1);
            rho_y1 = getRho(state_y1, nstate_y1);
            checkRange(testCase, rho_x1, rho_y1, a, b);
            assertNotEqual(testCase, rho_x1, rho_y1);

            c1(i, j) = ismembertol(nstate_x1.tree(r1).time, ...
                                   nstate_y1.tree(r1).time);
            o1(i, j) = getOverlap(state_x1, state_y1, a, b);

            % Already coupled
            [nstate_x2, nstate_y2, ~, ~, ~, ~, ~, ~] = RscaleTopTreeCoupled(...
                state_x2, state_y2, prior, del, deldel);

            rho_x2 = getRho(state_x2, nstate_x2);
            rho_y2 = getRho(state_y2, nstate_y2);
            checkRange(testCase, rho_x2, rho_y2, a, b);
            assertEqual(testCase, rho_x2, rho_y2);

            c2(i, j) = ismembertol(nstate_x2.tree(r2).time, ...
                                   nstate_y2.tree(r2).time);
            o2(i, j) = getOverlap(state_x2, state_y2, a, b);
        end
    end

    % Proportion of matching samples
    fprintf('Proportion of matching samples after %g trials\n', n_i * n_j);
    fprintf('Housekept trees\n');
    fprintf('  Theoretical          = %g\n', mean(o1(:)));
    fprintf('  RscaleTopTreeCoupled = %g\n', mean(c1(:)));
    fprintf('Coupled trees\n');
    fprintf('  Theoretical          = %g\n', mean(o2(:)));
    fprintf('  RscaleTopTreeCoupled = %g\n', mean(c2(:)));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function marginalTest(testCase)
    load('+RscaleTopTreeCoupledTest/cladePrior.mat', 'prior');
    rangeL = repelem(9:11, 3);
    n_i = length(rangeL);
    n_j = 2e3;
    del = 0.5;
    deldel = 1.5;
    a = del;
    b = del + deldel;
    [xc, yc, xm, ym] = deal(nan(n_i, n_j));
    for i = 1:n_i
        L = rangeL(i);
        [state_x, state_y] = RscaleTopTreeCoupledTest.housekeptStates(L);
        for j = 1:n_j
            % Coupled
            [nstate_xc, nstate_yc, ~, ~, ~, ~, ~, ~] = RscaleTopTreeCoupled(...
                state_x, state_y, prior, del, deldel);

            rho_xc = getRho(state_x, nstate_xc);
            rho_yc = getRho(state_y, nstate_yc);
            checkRange(testCase, rho_xc, rho_yc, a, b);

            xc(i, j) = rho_xc;
            yc(i, j) = rho_yc;

            % Marginal
            [nstate_xm, ~, ~, ~, ~] = RscaleTopTree(state_x, prior, del, deldel);
            [nstate_ym, ~, ~, ~, ~] = RscaleTopTree(state_y, prior, del, deldel);

            rho_xm = getRho(state_x, nstate_xm);
            rho_ym = getRho(state_y, nstate_ym);
            checkRange(testCase, rho_xm, rho_ym, a, b);

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
function ol = getOverlap(state_x, state_y, a, b)
    global LEAF;
    r = state_x.root;
    nu_x = progeny(state_x.tree, r, LEAF);
    nu_y = progeny(state_y.tree, r, LEAF);

    x_c = state_x.tree(r).time;
    y_c = state_y.tree(r).time;

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

function rho = getRho(state, nstate)
    r = state.root;
    t0 = min([state.tree(state.leaves).time]);
    rho = (nstate.tree(r).time - t0) / (state.tree(r).time - t0);
end

function checkRange(testCase, var_x, var_y, a, b)
    assertGreaterThanOrEqual(testCase, [var_x, var_y], a);
    assertLessThanOrEqual(testCase, [var_x, var_y], b);
end

% Setup and teardown functions
function setupOnce(testCase)
    RscaleTopTreeCoupledTest.setupOnce(testCase);
end

function tearDownOnce(testCase)
    RscaleTopTreeCoupledTest.tearDownOnce(testCase);
end

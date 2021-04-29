function tests = RscaleCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    n_i = 1e2;
    n_j = 1e2;
    a = 0.5;
    b = 2;
    [o1, o2] = deal(nan(n_i, 1));
    [c1, c2] = deal(nan(n_i, n_j));
    for i = 1:n_i
        [state_x1, state_y1] = housekeptStates();
        [state_x2, state_y2] = coupledStates();
        o1(i) = getOverlap(state_x1, state_y1, a, b);
        o2(i) = getOverlap(state_x2, state_y2, a, b);
        for j = 1:n_j
            [var_x1, var_y1] = RscaleCoupled(state_x1, state_y1, a, b);
            checkRange(testCase, var_x1, var_y1, a, b);
            assertNotEqual(testCase, var_x1, var_y1);

            n_x1 = getNewAge(var_x1, state_x1);
            n_y1 = getNewAge(var_y1, state_y1);
            c1(i, j) = ismembertol(n_x1, n_y1);

            [var_x2, var_y2] = RscaleCoupled(state_x2, state_y2, a, b);
            checkRange(testCase, var_x2, var_y2, a, b);
            assertEqual(testCase, var_x2, var_y2);

            n_x2 = getNewAge(var_x2, state_x2);
            n_y2 = getNewAge(var_y2, state_y2);
            c2(i, j) = ismembertol(n_x2, n_y2);
        end
    end

    % Proportion of matching samples
    fprintf('Proportion of matching samples after %g trials\n', n_i * n_j);
    fprintf('Housekept trees\n');
    fprintf('  theoretical                   = %g\n', mean(o1(:)));
    fprintf('  maximalCouplingUniformScaling = %g\n', mean(c1(:)));
    fprintf('Coupled trees\n');
    fprintf('  theoretical                   = %g\n', mean(o2(:)));
    fprintf('  maximalCouplingUniformScaling = %g\n', mean(c2(:)));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function marginalTest(testCase)
    n_i = 1e1;
    n_j = 1e4;
    a = 0.5;
    b = 2;
    [x0, y0] = deal(nan(n_i, 1));
    [xc, yc, xm, ym] = deal(nan(n_i, n_j));
    for i = 1:n_i
        [state_x, state_y] = housekeptStates();
        x0(i) = state_x.tree(state_x.root).time;
        y0(i) = state_y.tree(state_y.root).time;
        for j = 1:n_j
            [var_xc, var_yc] = RscaleCoupled(state_x, state_y, a, b);
            xc(i, j) = getNewAge(var_xc, state_x);
            yc(i, j) = getNewAge(var_yc, state_y);

            [nstate_x, ~, ~, ~, ~] = Rscale(state_x, a + rand * (b - a));
            [nstate_y, ~, ~, ~, ~] = Rscale(state_y, a + rand * (b - a));
            xm(i, j) = nstate_x.tree(nstate_x.root).time;
            ym(i, j) = nstate_y.tree(nstate_y.root).time;
        end
    end

    % Empirical CDFs
    [ecdf_xc, x_xc] = ecdf(xc(:));
    [ecdf_yc, x_yc] = ecdf(yc(:));

    [ecdf_xm, x_xm] = ecdf(xm(:));
    [ecdf_ym, x_ym] = ecdf(ym(:));

    % Actual CDFs
    f_cdf = @(t, z) mean((min(max(t, z * a), z * b) - z * a) ./ (z * (b - a)));
    cdf_xc = arrayfun(@(t) f_cdf(t, x0), x_xc(:));
    cdf_yc = arrayfun(@(t) f_cdf(t, y0), x_yc(:));
    cdf_xm = arrayfun(@(t) f_cdf(t, x0), x_xm(:));
    cdf_ym = arrayfun(@(t) f_cdf(t, y0), x_ym(:));

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
    
    fprintf('Difference from CDF for coupled and marginal samples');
    fprintf('%g draws for each of %g pairs of trees\n', n_j, n_i);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

%%%%

function [state_x, state_y] = housekeptStates()
    [L, theta] = getPars();
    state_x = dummyState(ExpTree(L, theta));
    state_y = housekeeping(state_x, dummyState(ExpTree(L, theta)));
end

function [state_x, state_y] = coupledStates()
    [L, theta] = getPars();
    [state_x, state_y] = deal(dummyState(ExpTree(L, theta)));
end

function state = dummyState(s)
    global ROOT
    % Same as makestate but without data and associated calculations
    % Assume catastrophes are already on tree
    state = tree2state(s);
    state.claderoot = [];
    state.cat = cellfun('length', {state.tree.catloc});
    state.cat = state.cat(:);
    state.ncat = sum(state.cat);
    state.length = TreeLength(state.tree, state.root);
    state.kappa = rand;
    state.root = find([state.tree.type] == ROOT);
end

% Helper functions
function ol = getOverlap(state_x, state_y, a, b)
    i = state_x.root;
    x = state_x.tree(i).time;
    y = state_y.tree(i).time;
    min_xy = min(x, y);
    max_xy = max(x, y);
    if b * min_xy < a * max_xy
        ol = 0;
    else
        ol = (b * min_xy - a * max_xy) / (max_xy * (b - a));
    end
end

function newage = getNewAge(var, state)
    newage = var * state.tree(state.root).time;
end

function checkRange(testCase, var_x, var_y, a, b)
    assertGreaterThanOrEqual(testCase, [var_x, var_y], a);
    assertLessThanOrEqual(testCase, [var_x, var_y], b);
end

function [L, theta] = getPars()
    L = 10;
    theta = 0.01;
end

% Setup and teardown functions
function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    global BORROWING MCMCCAT;
    testCase.TestData.BORROWING = BORROWING;
    testCase.TestData.MCMCCAT = MCMCCAT;
    BORROWING = 0;
    MCMCCAT = 0;
    warning('Assume no borrowing or catastrophes');
end

function teardownOnce(testCase)
    global BORROWING MCMCCAT
    BORROWING = testCase.TestData.BORROWING;
    MCMCCAT = testCase.TestData.MCMCCAT;
end

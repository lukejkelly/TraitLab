function tests = maximalCouplingUniformScalingTest
    % Unit-testing various maximal coupling functions
    tests = functiontests(localfunctions);
end

function directTest(testCase)
    xc = [0.9, 10, 1.1];
    yc = [0.5, 20, 1.11];
    for xy = [xc; yc]
        [indep, mc, mcl, mcus, ranges, ol] = directScaling(xy(1), xy(2));
        compareDistributions(testCase, indep, mc, mcl, mcus, ranges, ol);
    end
end

function indirectTest(testCase)
    xc = [0.9, 0.9];
    yc = [0.5, 0.1];
    for xy = [xc; yc]
        [indep, mc, mcl, mcus, ranges, ol] = indirectScaling(xy(1), xy(2));
        compareDistributions(testCase, indep, mc, mcl, mcus, ranges, ol);
    end
end

% Sampling functions
function [indep, mc, mcl, mcus, ranges, ol] = directScaling(xc, yc)
    % Sampling x = xc * U(a, b)
    if any([xc, yc] <= 0)
        error('Arguments must be positive');
    end

    [a, b, n, nn, rp, dp, ldp] = getParameters();

    % Independent samples
    indep = a * [xc, yc] + (b - a) * rand(n, 2) .* [xc, yc];

    % Maximal coupling using while loop and densities on original or log scale
    rx = @() rp(xc);
    ry = @() rp(yc);

    dx = @(x) dp(x, xc);
    dy = @(y) dp(y, yc);

    ldx = @(x) ldp(x, xc);
    ldy = @(y) ldp(y, yc);

    [mc_x, mc_y] = arrayfun(@(~) maximalCoupling(rx, dx, ry, dy), nn);
    mc = [mc_x, mc_y];

    [mcl_x, mcl_y] = arrayfun(@(~) maximalCouplingLog(rx, ldx, ry, ldy), nn);
    mcl = [mcl_x, mcl_y];

    % Maximal coupling without while loop
    [mcus_x, mcus_y] = arrayfun(@(~) maximalCouplingUniformScaling(xc, yc, ...
        a, b), nn);
    mcus = [mcus_x, mcus_y];

    % Endpoints of exact distribution
    ranges = [a, b]' * [xc, yc];

    % Overlap
    if (b * xc < a * yc || b * yc < a * xc)
        ol = 0;
    else
        if xc <= yc
            ol = (b * xc - a * yc) / (yc * (b - a));
        else
            ol = (b * yc - a * xc) / (xc * (b - a));
        end
    end
end

function [indep, mc, mcl, mcus, ranges, ol] = indirectScaling(xc, yc)
    % Sampling x = 1 - (1 - xc) * U(a, b)
    if any([xc, yc] <= 0 | [xc, yc] >= 1)
        error('Arguments must be between 0 and 1');
    end

    [a, b, n, nn, rp, dp, ldp] = getParameters();

    % Independent samples
    indep = 1 - b * (1 - [xc, yc]) + rand(n, 2) .* (1 - [xc, yc]) * (b - a);

    % Maximal coupling using while loop and densities on original or log scale
    rx = @() 1 - rp(1 - xc);
    ry = @() 1 - rp(1 - yc);

    dx = @(x) dp(1 - x, 1 - xc);
    dy = @(y) dp(1 - y, 1 - yc);

    ldx = @(x) ldp(1 - x, 1 - xc);
    ldy = @(y) ldp(1 - y, 1 - yc);

    [mc_x, mc_y] = arrayfun(@(~) maximalCoupling(rx, dx, ry, dy), nn);
    mc = [mc_x, mc_y];

    [mcl_x, mcl_y] = arrayfun(@(~) maximalCouplingLog(rx, ldx, ry, ldy), nn);
    mcl = [mcl_x, mcl_y];

    % Maximal coupling without while loop
    [mcus_x, mcus_y] = arrayfun(@(~) maximalCouplingUniformScaling(1 - xc, ...
        1 - yc, a, b), nn);
    mcus =  1 - [mcus_x, mcus_y];

    % Endpoints of exact distribution
    ranges = 1 - [b, a]' * (1 - [xc, yc]);

    % Overlap
    if (b * (1 - xc) < a * (1 - yc) || b * (1 - yc) < a * (1 - xc))
        ol = 0;
    else
        if xc <= yc
            ol = (b * (1 - yc) - a * (1 - xc)) / ((1 - xc) * (b - a));
        else
            ol = (b * (1 - xc) - a * (1 - yc)) / ((1 - yc) * (b - a));
        end
    end
end

% Comparison functions
function compareDistributions(testCase, indep, mc, mcl, mcus, ranges, ol)

    % Plotting empirical and actual CDFs
    clf();
    for j = 1:2
        subplot(1, 2, j);
        [n1, e1] = histcounts(indep(:, j), 101, 'Normalization', 'cdf');
        [n2, e2] = histcounts(mc(:, j), 101, 'Normalization', 'cdf');
        [n3, e3] = histcounts(mcl(:, j), 101, 'Normalization', 'cdf');
        [n4, e4] = histcounts(mcus(:, j), 101, 'Normalization', 'cdf');

        plot([e1; e2; e3; e4]', [zeros(1, 4); [n1; n2; n3; n4]'], ':', ...
            'LineWidth', 2); hold on;
        plot(ranges(:, j), [0, 1], '--', 'LineWidth', 2); hold off
        title(sprintf('U(%g, %g)', ranges(1, j), ranges(2, j)));
    end
    legend('indep', 'mc', 'mcl', 'mcus', 'exact');
    sgtitle('Exact and empirical CDFs');
    fmt = ' %-6.4g %-6.4g %-6.4g %-6.4g\n';
    fprintf('Sample ranges                 min_x  max_x  min_y  max_y\n');
    fprintf(['theoretical                  ', fmt], ranges);
    fprintf(['independent                  ', fmt], getRanges(indep));
    fprintf(['maximalCoupling              ', fmt], getRanges(mc));
    fprintf(['maximalCouplingLog           ', fmt], getRanges(mcl));
    fprintf(['maximalCouplingUniformScaling', fmt], getRanges(mcus));

    v1 = input('Do these CDFs and ranges match? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);

    % Proportion of matching samples
    fprintf('Proportion of matching samples\n');
    fprintf('theoretical                   = %g\n', ol);
    fprintf('maximalCoupling               = %g\n', propMatch(mc));
    fprintf('maximalCouplingLog            = %g\n', propMatch(mcl));
    fprintf('maximalCouplingUniformScaling = %g\n', propMatch(mcus));

    v2 = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v2 == 1);
end

% Helper functions
function [a, b, n, nn, rp, dp, ldp] = getParameters()
    a = 0.5;
    b = 2;
    n = 1e5;
    nn = (1:n)';
    rp = @(v) rpProto(v, a, b);
    dp = @(u, v) dpProto(u, v, a, b);
    ldp = @(u, v) ldpProto(u, v, a, b);
end

function u = rpProto(v, a, b)
    u = v * (a + rand * (b - a));
end

function d = dpProto(u, v, a, b)
    if (v * a <= u && u <= v * b)
        d = 1 / (v * (b - a));
    else
        d = 0;
    end
end

function l = ldpProto(u, v, a, b)
    if (v * a <= u && u <= v * b)
        l = -(log(v) + log(b - a));
    else
        l = -Inf;
    end
end

function p = propMatch(x)
    p = mean(x(:, 1) == x(:, 2));
end

function r = getRanges(xy)
    min_xy = min(xy);
    max_xy = max(xy);
    r = [min_xy(1), max_xy(1), min_xy(2), max_xy(2)];
end

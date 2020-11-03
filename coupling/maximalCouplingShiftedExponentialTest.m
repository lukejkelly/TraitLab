function tests = maximalCouplingShiftedExponentialTest
    % Unit-testing various maximal coupling functions
    tests = functiontests(localfunctions);
end

function samplingTest(testCase)
    a_p = [1, 10, 1000];
    a_q = [10, 10, 100];
    theta = [0.1, 0.01, 0.001];
    for i = 1:length(theta)
        compareDistributions(testCase, a_p(i), a_q(i), theta(i));
    end
end


% Sampling functions
function compareDistributions(testCase, a_p, a_q, theta)

    if any([a_p, a_q, theta] <= 0)
        error('Arguments must be positive');
    end

    n = 1e5;
    nn = (1:n)';

    rp = @() rpProto(a_p, theta);
    rq = @() rpProto(a_q, theta);

    ldp = @(x) ldpProto(x, a_p, theta);
    ldq = @(y) ldpProto(y, a_q, theta);

    % Independent samples
    indep_x = arrayfun(@(~) rp(), nn);
    indep_y = arrayfun(@(~) rq(), nn);
    indep = [indep_x, indep_y];

    % Maximal coupling using while loop and original/log densities
    [mcl_x, mcl_y] = arrayfun(@(~) maximalCouplingLog(rp, ldp, rq, ldq), nn);
    mcl = [mcl_x, mcl_y];

    % Maximal coupling without while loop
    [mcse_x, mcse_y] = arrayfun(@(~) maximalCouplingShiftedExponential(a_p, ...
        a_q, theta), nn);
    mcse = [mcse_x, mcse_y];

    % Overlap
    ol = exp(-theta * abs(a_p - a_q));


    % Plotting empirical and actual CDFs
    clf();
    for j = 1:2
        subplot(1, 2, j);
        [n1, e1] = histcounts(indep(:, j), 101, 'Normalization', 'cdf');
        [n2, e2] = histcounts(mcl(:, j), 101, 'Normalization', 'cdf');
        [n3, e3] = histcounts(mcse(:, j), 101, 'Normalization', 'cdf');

        semilogx([e1; e2; e3]', [zeros(1, 3); [n1; n2; n3]'], ':', ...
            'LineWidth', 2);
        if j == 1
            r = a_p;
        else
            r = a_q;
        end
        title(sprintf('ExpS(%g, %g)', theta, r));
    end

    legend('indep', 'mcl', 'mcse');
    suptitle(sprintf('Exact and empirical CDFs'));
    fmt = ' %-6.4g %-6.4g\n';
    fprintf('Sample minima:                min_x  min_y\n');
    fprintf(['actual                       ', fmt], a_p, a_q);
    fprintf(['independent                  ', fmt], min(indep_x), min(indep_y));
    fprintf(['maximalCouplingLog           ', fmt], min(mcl_x), min(mcl_y));
    fprintf(['maximalCouplingUniformScaling', fmt], min(mcse_x), min(mcse_y));

    v1 = input('Do these CDFs and ranges match? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);

    % Proportion of matching samples
    fprintf('Proportion of matching samples\n');
    fprintf('theoretical                       = %g\n', ol);
    fprintf('maximalCouplingLog                = %g\n', propMatch(mcl));
    fprintf('maximalCouplingShiftedExponential = %g\n', propMatch(mcse));

    v2 = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v2 == 1);
end

% Helper functions
function u = rpProto(a, theta)
    u = a - log(rand) / theta;
end

function l = ldpProto(u, a, theta)
    if (u >= a)
        l = log(theta) - theta * (u - a);
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

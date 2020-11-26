function tests = maximalCouplingUniformShiftedExponentialTest
    % Unit-testing various maximal coupling functions
    tests = functiontests(localfunctions);
end

function noOverlapTest(testCase)
    compareDistributions(testCase, 1, 2, 2, 1);
end

function partialOverlap1Test(testCase)
    compareDistributions(testCase, 1.5, 3, 2, 1);
end

function partialOverlap2Test(testCase)
    compareDistributions(testCase, 2.5, 3, 2, 1);
end

function partialOverlap3Test(testCase)
    compareDistributions(testCase, 3, 5, 2, 1);
end

function compareDistributions(testCase, a, b, c, theta)

    n = 5e5;
    [x, y] = deal(nan(n, 1));
    for i = 1:n
        [x(i), y(i)] = maximalCouplingUniformShiftedExponential(a, b, c, theta);
    end

    nBins = 1e2;
    subplot(1, 2, 1);
    histogram(x, nBins, 'Normalization', 'pdf');
    hold on;
    plot([a, b], repmat(1 / (b - a), 1, 2), 'LineWidth', 2);
    hold off;
    legend('Obs', 'Exp');
    title(sprintf('x: uniform(%g, %g)', a, b));

    subplot(1, 2, 2)
    h = histogram(y, nBins, 'Normalization', 'pdf');
    hold on;
    plot(h.BinEdges, theta * exp(-theta * (h.BinEdges - c)), 'LineWidth', 2);
    hold off;
    legend('Obs', 'Exp');
    title(sprintf('y: exponential(%g) shifted %g', theta, c));

    oObs = mean(x == y);
    if b <= c
        oExp = 0;
    else
        % min(p, q) = p for x <= z, q otherwise; intersection possibly outside
        % range but doesn't matter here
        z = c + log(theta * (b - a)) / theta;
        if a <= c
            % supp min(p, q) = [c, b]
            oExp = max(0, z - c) / (b - a) + exp(-theta * (max(z - c, 0))) ...
                - exp(-theta * (b - c));
        else
            % c < a < b
            % supp min(p, q) = [a, b]
            oExp = max(0, z - a) / (b - a) + exp(-theta * (max(z, a) - c)) ...
                - exp(-theta * (b - c));
        end
    end
    fprintf('Proportion of coupling events\n');
    fprintf('Observed:   %g\nExpected:   %g\nDifference: %g\n', ...
             oObs, oExp, oObs - oExp);

    v = input('Do the CDFs and coupling proportions match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end

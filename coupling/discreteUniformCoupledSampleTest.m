function tests = discreteUniformCoupledSampleTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    N = 10;
    n_i = 5;
    n_j = 2e4;

    fprintf('Proportion of matching samples in each of %g trials\n', n_j);
    fprintf('%-4s%-8s%-8s%-8s\n', 'i', 'Obs', 'Exp1', 'Exp2');

    for i = 1:n_i

        n_x = 1 + floor(rand * N);
        s_x = randsample(N, n_x);

        n_y = 1 + floor(rand * N);
        s_y = randsample(N, n_y);

        [x, y] = deal(nan(1, n_j));
        for j = 1:n_j
            [x(j), y(j)] = discreteUniformCoupledSample(s_x, s_y);
        end

        pObs = cumsum(mean(x(:) == 1:N));
        pExp = cumsum(arrayfun(@(x) exp(discreteUniformLogProb(x, s_x)), 1:N));
        pAct = cumsum(ismember(1:N, s_x) / n_x);

        qObs = cumsum(mean(y(:) == 1:N));
        qExp = cumsum(arrayfun(@(y) exp(discreteUniformLogProb(y, s_y)), 1:N));
        qAct = cumsum(ismember(1:N, s_y) / n_y);

        subplot(n_i, 2, 2 * i - 1);
        plot(1:N, pObs, 'o', 1:N, pExp, 'x', 1:N, pAct, '^');
        xlabel('x');
        ylabel('CDF');
        title(sprintf('n = %d', n_x));

        subplot(n_i, 2, 2 * i);
        plot(1:N, qObs, 'o', 1:N, qExp, 'x', 1:N, qAct, '^');
        xlabel('y');
        ylabel('CDF');
        title(sprintf('n = %d', n_y));

        if i == n_i
            legend('ECDF', 'CDF1', 'CDF2', 'Location', 'southeast', ...
                   'Orientation', 'horizontal');
        end


        oObs = mean(x == y);
        oExp = sum(exp(arrayfun(...
            @(z) min(cellfun(@(s) discreteUniformLogProb(z, s), ...
                             {s_x, s_y})), ...
            1:N)));
        oAct = length(intersect(s_x, s_y)) / max(n_x, n_y);
        fprintf('%-4d%-8.4g%-8.4g%-8.4g\n', i, oObs, oExp, oAct);
    end
    fprintf('\nFigure displays marginal distributions\n');
    v = input('Do these distributions and proportions match? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function setupOnce(~)
    clf;
end

function teardownOnce(~)
    close;
end

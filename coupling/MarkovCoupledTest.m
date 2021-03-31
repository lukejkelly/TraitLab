function tests = MarkovCoupledTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    % Global variables from pausing runmcmcCoupled and entering an otherwise
    % empty workspace and saving what remained
    load('coupling/MarkovTestGlobal-20210330.mat');
    testCase.TestData.global_vars = whos('global');
end

function teardownOnce(testCase)
    cellfun(@clear, {testCase.TestData.global_vars.name});
end

function coupledTest(testCase)
    % States from runmcmcCoupled main loop just after successful coupling
    load('coupling/MarkovTestCoupled-20210330.mat', 'mcmc', 'model', ...
         'state_x', 'state_y');
    mcmc.subsample = 1;
    assertTrue(testCase, checkCoupling(state_x, state_y));
    n = 1e3;
    fprintf('\n');
    for i = 1:n
        if mod(i, floor(n / 10)) == 0
            fprintf('%.02g ', i / n);
        end
        [state_x, state_y, pa_x, pa_y] = MarkovCoupled(mcmc, model, ...
            state_x, state_y, 0);
        assertTrue(testCase, checkCoupling(state_x, state_y));
        assertEqual(testCase, pa_x, pa_y);
    end
end

function uncoupledTest(testCase)
    % Already checked during simulations
    % States from start of runmcmcCoupled main loop before coupling
    load('coupling/MarkovTestUncoupled-20210331.mat', 'mcmc', 'model', ...
         'state_x', 'state_y');
    assertFalse(testCase, checkCoupling(state_x, state_y));
    % Uneven move distribution
    [state_x, ~] = Markov(mcmc, model, state_x);
    [state_y, ~] = Markov(mcmc, model, state_y);
    compareMarginalCoupled(testCase, mcmc, model, state_x, state_y);
    % Even move distribution
    [state_x, ~] = Markov(mcmc, model, state_x);
    [state_y, ~] = Markov(mcmc, model, state_y);
    mcmc.update.move(find(mcmc.update.move)) = 1 / sum(mcmc.update.move > 0);
    mcmc.update.cmove = cumsum(mcmc.update.move);
    compareMarginalCoupled(testCase, mcmc, model, state_x, state_y);
end

function compareMarginalCoupled(testCase, mcmc, model, state_x, state_y)
    % Already checked during simulations
    % States from start of runmcmcCoupled main loop before coupling
    mcmc.subsample = 1;
    n = 1e4;
    [t_x1, t_y1, t_x2, t_y2] = deal(nan(n, 3));
    [a_x1, a_y1, a_x2, a_y2, b_x1, b_y1, b_x2, b_y2] = deal(nan(n, 1));
    fprintf('\n');
    for i = 1:n
        if mod(i, floor(n / 10)) == 0
            fprintf('%.02g ', i / n);
        end
        % Run coupled and marginal kernels
        [state_x1, state_y1, pa_x1, pa_y1] = MarkovCoupled(mcmc, model, ...
            state_x, state_y);
        [state_x2, pa_x2] = Markov(mcmc, model, state_x);
        [state_y2, pa_y2] = Markov(mcmc, model, state_y);
        % Store prior, likelihood and root time
        t_x1(i, :) = stateOutputs(state_x1);
        t_y1(i, :) = stateOutputs(state_y1);
        t_x2(i, :) = stateOutputs(state_x2);
        t_y2(i, :) = stateOutputs(state_y2);
        % % Store move attempts and successes
        [a_x1(i), b_x1(i)] = moveOutputs(pa_x1);
        [a_y1(i), b_y1(i)] = moveOutputs(pa_y1);
        [a_x2(i), b_x2(i)] = moveOutputs(pa_x2);
        [a_y2(i), b_y2(i)] = moveOutputs(pa_y2);
    end
    d_x1 = logical(b_x1);
    d_y1 = logical(b_y1);
    d_x2 = logical(b_x2);
    d_y2 = logical(b_y2);
    % Comparing sample distributions
    labs = {'log prior', 'log likelihood', 'root time'};
    clf;
    for j = 1:3
        subplot(2, 3, j);
        plotStateOutputs(t_x1(d_x1, j), t_x2(d_x2, j));
        title(sprintf('x : %s', labs{j}));
        subplot(2, 3, 3 + j);
        plotStateOutputs(t_y1(d_y1, j), t_y2(d_y2, j));
        title(sprintf('y : %s', labs{j}));
    end
    v1 = input('Do the CDFs for accepted moves match? Reply 1 for yes... ');
    assertEqual(testCase, v1, 1);
    % Compare counts of attempted moves and success rates
    clf;
    subplot(2, 2, 1);
    plotMoveCounts(a_x1, b_x1, a_x2, b_x2);
    title('x move counts');
    subplot(2, 2, 2);
    plotMoveProportions(a_x1, b_x1, a_x2, b_x2);
    title('x move successes');
    subplot(2, 2, 3);
    plotMoveCounts(a_y1, b_y1, a_y2, b_y2);
    title('y move counts');
    subplot(2, 2, 4);
    plotMoveProportions(a_y1, b_y1, a_y2, b_y2);
    title('y move successes');
    v2 = input('Do the move statistics match? Reply 1 for yes... ');
    assertEqual(testCase, v2, 1);
    % Posteriors for successful moves by type
    clf;
    m = unique([a_x1(d_x1); a_y1(d_y1); a_x2(d_x2); a_y2(d_y2)]);
    for i = 1:length(m)
        subplot(length(m), 2, 2 * i - 1);
        plotPosteriorByMove(t_x1, a_x1, d_x1, t_x2, a_x2, d_x2, m(i));
        if i == 1
            title('x posterior');
        end
        subplot(length(m), 2, 2 * i);
        plotPosteriorByMove(t_y1, a_y1, d_y1, t_y2, a_y2, d_y2, m(i));
        if i == 1
            title('y posterior');
        end
    end
    v3 = input('Do the posterior CDFs by move match? Reply 1 for yes... ');
    assertEqual(testCase, v3, 1);
end

function t = stateOutputs(state)
    t = [state.logprior, state.loglkd, state.tree(state.root).time];
end

function [a, b] = moveOutputs(pa)
    a = find(~isnan(pa));
    b = pa(a);
end

function plotStateOutputs(t1, t2)
    [~, edges] = histcounts([t1(:); t2(:)], 50);
    n1 = histcounts(t1, edges, 'Normalization', 'cdf');
    n2 = histcounts(t2, edges, 'Normalization', 'cdf');
    plot(edges, [0, n1], 'b-', 'LineWidth', 2); hold on;
    scatter(edges, [0, n2], '*r', 'MarkerFaceAlpha', 0.5); hold off;
    legend({sprintf('C %i', length(t1)), sprintf('M %i', length(t2))}, ...
           'Location', 'SouthEast');
    ylabel('ecdf');
    axis('tight');
end

function plotMoveCounts(a1, b1, a2, b2)
    m = unique([a1, a2]);
    n_a1 = arrayfun(@(i) sum(a1 == i), m);
    n_b1 = arrayfun(@(i) sum(b1(a1 == i)), m);
    n_a2 = arrayfun(@(i) sum(a2 == i), m);
    n_b2 = arrayfun(@(i) sum(b2(a2 == i)), m);
    x = categorical(m);
    scatter(x, n_a1, 'or', 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    scatter(x, n_b1, 'sr', 'filled', 'MarkerFaceAlpha', 0.5);
    scatter(x, n_a2, '*b', 'MarkerFaceAlpha', 0.5);
    scatter(x, n_b2, '^b', 'MarkerFaceAlpha', 0.5);
    hold off;
    legend({'coupled attempted', 'coupled succeeded', ...
            'marginal attempted', 'marginal succeeded'}, ...
            'Location', 'southoutside', 'Orientation', 'horizontal');
    axis('tight');
    xlabel('move');
    ylabel('count')
end

function plotMoveProportions(a1, b1, a2, b2)
    m = unique([a1, a2]);
    n_a1 = arrayfun(@(i) sum(a1 == i), m);
    n_b1 = arrayfun(@(i) sum(b1(a1 == i)), m);
    n1 = n_b1 ./ n_a1;
    n_a2 = arrayfun(@(i) sum(a2 == i), m);
    n_b2 = arrayfun(@(i) sum(b2(a2 == i)), m);
    n2 = n_b2 ./ n_a2;
    x = categorical(m);
    scatter(x, n1, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    scatter(x, n2, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
    hold off;
    legend({'coupled', 'marginal'}, 'Location', 'southoutside', ...
           'Orientation', 'horizontal');
    axis('tight');
    xlabel('move');
    ylabel('success rate')
end

function plotPosteriorByMove(t1, a1, d1, t2, a2, d2, mi)
    plotStateOutputs(sum(t1(a1 == mi & d1, 1:2), 2), ...
                     sum(t2(a2 == mi & d2, 1:2), 2));
    xlabel(moveNames(mi));
end

function mn = moveNames(i)
    mns = {'slide', 'exchange (local)', 'exchange (wide)', 'spr (local)', ...
           'spr (wide)', 'scale (tree)', 'scale (subtree)', 'mu', [], [], ...
           'vary leaves', 'scale (top)', 'catastrophe (add)', ...
           'catastrophe (delete)', 'rho', 'kappa', [], 'catastrophe (move)', ...
           'xi (one)', 'xi (all)', 'beta'};
    mn = mns{i};
end

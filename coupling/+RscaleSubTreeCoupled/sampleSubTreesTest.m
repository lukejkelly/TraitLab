function tests = sampleSubTreesTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global LEAF;
    rangeL = 9:12;
    nSamp = 1e4;
    for i = 1:length(rangeL)
        L = rangeL(i);
        [state_x, state_y] ...
            = RscaleSubTreeCoupled.unitTests.coupledStates(L, 1);
        % Expected output
        nu = progeny(state_x.tree, state_x.root, LEAF);
        pExp = nu(2, :) / sum(nu(2, :));
        % Observed frequencies
        [i_x, i_y] = arrayfun(...
            @(~) RscaleSubTreeCoupled.sampleSubTrees(state_x, state_y), ...
            1:nSamp);
        assertTrue(testCase, all(i_x == i_y));
        pObs = arrayfun(@(z) mean(i_x == z), nu(1, :));
        subplot(2, 2, i);
        scatter(pObs, pExp, 'MarkerEdgeAlpha', 0.75);
        refline(1, 0);
        xlabel('Expected proportion');
        ylabel('Observed proportion');
        title(sprintf('L = %d leaves', L));
    end
    fprintf('Observed and expected proportions from %g samples\n', nSamp);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function marginalTest(testCase)
    global LEAF;
    rangeL = 9:12;
    nSamp = 1e4;
    for i = 1:length(rangeL)
        L = rangeL(i);
        [state_x, state_y] ...
            = RscaleSubTreeCoupled.unitTests.housekeptStates(L, 1);
        % Expected output
        nu_x = sortrows(progeny(state_x.tree, state_x.root, LEAF)')';
        nu_y = sortrows(progeny(state_y.tree, state_y.root, LEAF)')';
        xExp = nu_x(2, :) / sum(nu_x(2, :));
        yExp = nu_y(2, :) / sum(nu_y(2, :));
        % Observed output
        [i_x, i_y, nu_xObs, nu_yObs] = arrayfun(...
            @(~) RscaleSubTreeCoupled.sampleSubTrees(state_x, state_y), ...
            1:nSamp, 'UniformOutput', false);
        xObs = arrayfun(@(z) mean([i_x{:}] == z), nu_x(1, :));
        yObs = arrayfun(@(z) mean([i_y{:}] == z), nu_y(1, :));
        % Compare marginal frequencies
        subplot(4, 2, 2 * i - 1);
        scatter(xObs, xExp, 'MarkerEdgeAlpha', 0.75); hold on;
        scatter(yObs, yExp, 'MarkerEdgeAlpha', 0.75); hold off;
        refline(1, 0);
        xlabel('Expected');
        ylabel('Observed');
        legend('x', 'y', 'Location', 'SouthEast');
        title(sprintf('L = %d : marginal', L));
        % Compare coupling frequencies
        zExp = min(xExp, yExp);
        zObs = arrayfun(@(z) mean([i_x{:}] == z & [i_x{:}] == [i_y{:}]), ...
                        nu_x(1, :));
        subplot(4, 2, 2 * i);
        scatter(zObs, zExp, 'MarkerEdgeAlpha', 0.75); hold on;
        refline(1, 0);
        xlabel('Expected');
        ylabel('Observed');
        title(sprintf('L = %d : coupling', L));
        % Compare progeny
        cellfun(@(z_i, z_nu) assertTrue(testCase, ...
            all(progeny(state_x.tree, z_i, LEAF) == z_nu, 'all')), i_x, nu_xObs);
        cellfun(@(z_i, z_nu) assertTrue(testCase, ...
            all(progeny(state_y.tree, z_i, LEAF) == z_nu, 'all')), i_y, nu_yObs);
    end
    fprintf('Observed and expected proportions from %g samples\n', nSamp);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function setupOnce(testCase)
    RscaleSubTreeCoupled.unitTests.setupOnce(testCase);
end

function setup(~)
    clf;
end

function teardownOnce(testCase)
    RscaleSubTreeCoupled.unitTests.teardownOnce(testCase);
    close;
end

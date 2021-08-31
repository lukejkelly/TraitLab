function tests = rDiscreteTest
    tests = functiontests(localfunctions);
end

function overallTest(testCase)
    global LEAF;
    nSamp = 1e4;
    rangeL = 5:8;
    for i = 1:length(rangeL)
        L = rangeL(i);
        state = RscaleSubTreeCoupled.unitTests.dummyState(ExpTree(L, 1));
        subplot(2, 2, i);
        hold on;
        for j = state.nodes
            nu = progeny(state.tree, j, LEAF);
            x = arrayfun(@(~) RscaleSubTreeCoupled.rDiscrete(nu), 1:nSamp);
            pObs = arrayfun(@(z) mean(x == z), nu(1, :));
            pExp = exp(RscaleSubTreeCoupled.lpDiscrete(nu(1, :), nu));
            scatter(pExp, pObs, 'MarkerEdgeAlpha', 0.75);
        end
        hold off;
        refline(1, 0);
        xlabel('Expected proportion');
        ylabel('Observed proportion');
        title(sprintf('L = %d leaves', L));
    end
    fprintf('Observed and expected proportions from %g samples\n', nSamp);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function setupOnce(testCase)
    RscaleSubTreeCoupled.unitTests.setupOnce(testCase);
    figure;
end

function teardownOnce(testCase)
    RscaleSubTreeCoupled.unitTests.teardownOnce(testCase);
    close;
end

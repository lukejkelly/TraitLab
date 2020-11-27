function tests = BchooseCoupledMaximalTest
    % Compare marginal behavuour to core/Bchoose
    % Mostly a duplicate of BchooseCoupledMaximal.sampleMarginalTest
    tests = functiontests(localfunctions);
end

function cladesNoMoveNarrowTest(testCase)
    compareAllDistributions(testCase, 'No', 'Narrow', 'No');
    compareAllDistributions(testCase, 'No', 'Narrow', 'Yes');
end

function cladesYesMoveNarrowTest(testCase)
    compareAllDistributions(testCase, 'Yes', 'Narrow', 'No');
    compareAllDistributions(testCase, 'Yes', 'Narrow', 'Yes');
end

function cladesNoMoveWideTest(testCase)
    compareAllDistributions(testCase, 'No', 'Wide', 'No');
    compareAllDistributions(testCase, 'No', 'Wide', 'Yes');
end

function cladesYesMoveWideTest(testCase)
    compareAllDistributions(testCase, 'Yes', 'Wide', 'No');
    compareAllDistributions(testCase, 'Yes', 'Wide', 'Yes');
end

% Helper functions
function compareAllDistributions(testCase, clades, moveType, flip)
    % For now we look at overall sampling behaviour rather than comparing
    % distributions at each branch
    [state_x, state_y, mt, prior, nReps] = getParams(clades, moveType, flip);
    [xObs, yObs, xExp, yExp] = deal(nan(nReps, 5));

    for r = 1:nReps
        THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();

        % Coupled samples
        [iObs, jObs_x, jObs_y, kObs_x, kObs_y, newageObs_x, newageObs_y, ...
         logqObs_x, logqObs_y] ...
            = BchooseCoupledMaximal(state_x, state_y, mt, THETA, prior);

        xObs(r, :) = storeOutputs(iObs, jObs_x, kObs_x, newageObs_x, logqObs_x);
        yObs(r, :) = storeOutputs(iObs, jObs_y, kObs_y, newageObs_y, logqObs_y);

        % Marginal samples
        [iExp_x, jExp_x, kExp_x, newageExp_x, logqExp_x] ...
            = Bchoose(state_x, mt, THETA, prior);
        [iExp_y, jExp_y, kExp_y, newageExp_y, logqExp_y] ...
            = Bchoose(state_y, mt, THETA, prior);

        xExp(r, :) = storeOutputs(iExp_x, jExp_x, kExp_x, newageExp_x, logqExp_x);
        yExp(r, :) = storeOutputs(iExp_y, jExp_y, kExp_y, newageExp_y, logqExp_y);
    end

    parNames = {'i', 'j', 'k', 'newage', 'logq'};
    for i = 1:5
        subplot(5, 2, 2 * i - 1);
        drawECDFs(xObs(:, i), xExp(:, i), sprintf('x: %s', parNames{i}));
        subplot(5, 2, 2 * i);
        drawECDFs(yObs(:, i), yExp(:, i), sprintf('y: %s', parNames{i}));
    end
    fprintf('clades: %s\nmove type: %s\nflip: %s\nsamples: %g\n', ...
        clades, moveType, flip, nReps);
    fprintf('Histograms do not reach 1 if any moves fail\n');
    v = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end

function drawECDFs(xObs, xExp, figTitle)
    nBins = 50;
    [edges, nObs, nExp] = deal(zeros(nBins + 1, 1));
    [~, edges(:)] = histcounts([xObs, xExp], nBins);
    nObs(2:end) = histcounts(xObs, edges, 'Normalization', 'cdf');
    nExp(2:end) = histcounts(xExp, edges, 'Normalization', 'cdf');
    yyaxis left;
    plot(edges, nObs, '-b', edges, nExp, ':g', 'LineWidth', 2);
    yyaxis right;
    plot(edges, nObs - nExp, '-.r', 'LineWidth', 2);
    legend('Obs', 'Exp', 'Diff', 'Location', 'SouthEast');
    axis('tight');
    title(figTitle);
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function [state_x, state_y, mt, prior, nReps] = getParams(clades, moveType, flip)
    global NARROW WIDE
    [~, state_x] = BchooseCoupledMaximal.state10a(clades);
    [~, state_y] = BchooseCoupledMaximal.state10b(clades);
    if strcmp(flip, 'Yes')
        [state_x, state_y] = deal(state_y, state_x);
    end
    state_y = housekeeping(state_x, state_y);
    switch moveType
    case 'Narrow'
        mt = NARROW;
    case 'Wide'
        mt = WIDE;
    otherwise
        error('moveType must be ''NARROW'' or ''WIDE''');
    end
    prior.isclade = strcmp(clades, 'yes');
    nReps = 5e5;
end

function x = storeOutputs(i, j, k, newage, logq)
    if isempty(newage)
        [newage, logq] = deal(nan);
    end
    x = [i, j, k, newage, logq];
end

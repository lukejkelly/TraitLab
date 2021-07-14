function tests = BchooseCoupledTest
    % Compare marginal behavuour to core/Bchoose
    % Mostly a duplicate of BchooseCoupled.sampleMarginalTest
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
        THETA = BchooseCoupled.sampleCoupling.sampleTheta();

        % Coupled samples
        [iObs, jObs_x, jObs_y, kObs_x, kObs_y, newageObs_x, newageObs_y, ...
         logqObs_x, logqObs_y] = BchooseCoupled(state_x, state_y, mt, THETA, prior);

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
        drawECDFs(xObs(:, i), xExp(:, i), sprintf('x: %s', parNames{i}), i < 4);
        subplot(5, 2, 2 * i);
        drawECDFs(yObs(:, i), yExp(:, i), sprintf('y: %s', parNames{i}), i < 4);
    end
    fprintf('clades: %s\nmove type: %s\nflip: %s\nsamples: %g\n', ...
            clades, moveType, flip, nReps);
    fprintf('Histograms do not reach 1 if any moves fail\n');
    v1 = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v1, 1);

    drawIndividualECDFs(state_x, xObs, xExp, 'x');
    drawIndividualECDFs(state_y, yObs, yExp, 'y');
    v2 = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v2, 1);
end

function drawECDFs(xObs, xExp, figTitle, binIntegers)
    if binIntegers
        [~, edges] = histcounts([xObs, xExp], 'BinMethod', 'integers');
    else
        [~, edges] = histcounts([xObs, xExp], 50);
    end
    edges = edges(:);
    nObs = histcounts(xObs, edges, 'Normalization', 'cdf');
    nExp = histcounts(xExp, edges, 'Normalization', 'cdf');
    binCentres = conv(edges, [0.5, 0.5], 'valid');
    yyaxis left;
    plot(binCentres, nObs, '-b', binCentres, nExp, ':g', 'LineWidth', 2);
    ylim([0, 1]);
    yyaxis right;
    plot(binCentres, nObs - nExp, '-.r', 'LineWidth', 2);
    legend('Obs', 'Exp', 'Diff', 'Location', 'SouthEast');
    axis('tight');
    title(figTitle);
end

function drawIndividualECDFs(state, xObs, xExp, figVar)
    % Plot distribution of t_pa(i)' against pa(i)
    pObs = [state.tree(xObs(:, 1)).parent];
    tObs = xObs(:, 4);
    pExp = [state.tree(xExp(:, 1)).parent];
    tExp = xExp(:, 4);

    np = state.NS - 1;
    uObs = unique(pObs);
    uExp = unique(pExp);
    if any(np ~= [length(uObs), length(uExp)]) || ~all(uObs == uExp) ...
           || ~all(sort(state.nodes) == uObs)
        error('Coverage does not overlap');
    end

    if strcmp(figVar, 'x')
        figInds = 1:2:(2 * np - 1);
    else
        figInds = 2:2:(2 * np);
    end
    for i = 1:np
        subplot(np, 2, figInds(i));
        tObs_i = tObs(pObs == uObs(i));
        tExp_i = tExp(pExp == uExp(i));
        [~, edges] = histcounts([tObs_i; tExp_i], 50);
        nObs = histcounts(tObs_i, edges, 'Normalization', 'cdf');
        nExp = histcounts(tExp_i, edges, 'Normalization', 'cdf');
        binCentres = conv(edges, [0.5, 0.5], 'valid');
        yyaxis left;
        plot(binCentres, nObs, '-b', binCentres, nExp, ':g', 'LineWidth', 2);
        ylim([0, 1]);
        yyaxis right;
        plot(binCentres, nObs - nExp, '-.r', 'LineWidth', 2);
        legend('Obs', 'Exp', 'Diff', 'Location', 'SouthEast');
        axis('tight');
        title(sprintf('%s: pa(i) = %d', figVar, uObs(i)));
    end
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
    global BORROWING
    if isempty(BORROWING)
        BORROWING = 0;
    end
    global MCMCCAT
    if isempty(MCMCCAT)
        MCMCCAT = 0;
    elseif MCMCCAT == 1
        error('Code has not been checked for catastrophes');
    end
end

function [state_x, state_y, mt, prior, nReps] = getParams(clades, moveType, flip)
    global NARROW WIDE
    [~, state_x] = BchooseCoupled.state10a(clades);
    [~, state_y] = BchooseCoupled.state10b(clades);
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
    nReps = 1e4;
end

function x = storeOutputs(i, j, k, newage, logq)
    if isempty(newage)
        [newage, logq] = deal(nan);
    end
    x = [i, j, k, newage, logq];
end

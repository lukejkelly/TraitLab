function tests = BchooseCoupledTest
    % Compare marginal behaviour to core/Bchoose
    % Mostly a duplicate of Bchoose.sampleMarginalTest
    % Catastrophes are ignored except for the final move
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

function catastropheTest(testCase)
    % check catastrophe-related outputs only
    global BORROWING MCMCCAT NARROW WIDE
    M = 10;
    N = 10;
    for BORROWING = [0, 1]
        for MCMCCAT = [0, 1]
            for mt = [NARROW, WIDE]
                for m = 1:M
                    L = 5 + ceil(rand * 5);
                    [state, theta, prior] = dummyState(L);
                    for n = 1:N
                        [~, ~, ~, ~, ~, newage_x, newage_y, ~, ~, ...
                         ncat_x, ncat_y, cat_x, cat_y, loc_x, loc_y] ...
                            = BchooseCoupled(state, state, mt, theta, prior);

                        if isempty(newage_x) || ~MCMCCAT
                            assertEmpty(testCase, ncat_x);
                            assertEmpty(testCase, cat_x);
                            assertEmpty(testCase, loc_x);
                        else
                            assertGreaterThanOrEqual(testCase, ncat_x, 0);
                            assertNotEmpty(testCase, cat_x);
                            if BORROWING
                                assertNotEmpty(testCase, loc_x);
                            else
                                assertEmpty(testCase, loc_x);
                            end
                        end
                        if isempty(newage_y) || ~MCMCCAT
                            assertEmpty(testCase, ncat_y);
                            assertEmpty(testCase, cat_y);
                            assertEmpty(testCase, loc_y);
                        else
                            assertGreaterThanOrEqual(testCase, ncat_y, 0);
                            assertNotEmpty(testCase, cat_y);
                            if BORROWING
                                assertNotEmpty(testCase, loc_y);
                            else
                                assertEmpty(testCase, loc_y);
                            end
                        end
                        if MCMCCAT && ~isempty(newage_x) && ~isempty(newage_y)
                            assertEqual(testCase, ncat_x, ncat_y);
                            assertEqual(testCase, cat_x, cat_y);
                            if BORROWING
                                assertEqual(testCase, loc_x, loc_y);
                            end
                        end
                    end
                end
            end
        end
    end
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

    parNames = {'$ i $', '$ j $', '$ \mathrm{pa}(j) $', '$ t_i $', 'logq'};
    for i = 1:4 %5
        subplot(4, 2, 2 * i - 1);
        drawECDFs(xObs(:, i), xExp(:, i), sprintf('x: %s', parNames{i}), i < 4);
        subplot(4, 2, 2 * i);
        drawECDFs(yObs(:, i), yExp(:, i), sprintf('y: %s', parNames{i}), i < 4);
    end
    fprintf('\nclades: %s\nmove type: %s\nflip: %s\nsamples: %g\n', ...
            clades, moveType, flip, nReps);
    fprintf('Histograms do not reach 1 if any moves fail\n');
    v1 = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v1, 1);

    clf;
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
    ylabel('ECDF');
    ylim([0, 1]);
    yyaxis right;
    plot(binCentres, nObs - nExp, '-.r', 'LineWidth', 2);
    legend('Coupled', 'Marginal', 'Difference', 'Location', 'northoutside', ...
           'Orientation', 'horizontal');
    axis('tight');
    xlabel(figTitle, 'Interpreter', 'latex');
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
        ylabel('ECDF');
        ylim([0, 1]);
        yyaxis right;
        plot(binCentres, nObs - nExp, '-.r', 'LineWidth', 2);
        legend('Obs', 'Exp', 'Diff', 'Location', 'SouthEast');
        axis('tight');
        title(sprintf('%s: pa(i) = %d', figVar, uObs(i)));
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

function [state, theta, prior] = dummyState(L)
    global ROOT
    theta = rand * 1e-2;

    s = ExpTree(L, theta);
    clade = synthclades(s, ceil(rand * L / 2), 2, 1 - rand^3);
    rootmax = (1 + rand) * s([s.type] == ROOT).time;
    prior = unitTests.clade2prior(clade, rootmax);

    state = unitTests.dummyState(s);
    state = UpdateClades(state, [state.leaves, state.nodes], ...
                         size(prior.clade, 2));
    for c = 1:poissrnd(4)
        state = AddCat(state);
    end
    state.rho = state.ncat / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    global MCMCCAT VARYRHO BORROWING
    unitTests.setupOnce(testCase);
    testCase.TestData.VARYRHO = VARYRHO;
    MCMCCAT = 0;
    VARYRHO = 1;
    BORROWING = 0;
    clf;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

function teardown(~)
    close;
end

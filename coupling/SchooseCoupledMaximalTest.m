function tests = SchooseCoupledMaximalTest
    tests = functiontests(localfunctions);
end

function runTest(testCase)
    global BORROWING MCMCCAT
    for BORROWING = [0, 1]
        for MCMCCAT = [0, 1]
            parameterTests(testCase);
            remainCoupledTests(testCase);
            marginalTests(testCase);
            coupledTests(testCase);
        end
    end
end
function parameterTests(testCase)
    [state_x, state_y] = dummyStates();

    assertEqual(testCase, state_x.root, state_y.root);
    assertEqual(testCase, state_x.nodes, state_y.nodes);
    assertEqual(testCase, state_x.NS - 1, length(state_x.nodes));

    visited = zeros(1, state_x.NS - 1);

    while all(visited < 1e2)
        [i, newage_x, newage_y, logq_x, logq_y] ...
            = SchooseCoupledMaximal(state_x, state_y);
        visited(state_x.nodes == i) = visited(state_x.nodes == i) + 1;

        iT_x = state_x.tree(i).time;
        iT_y = state_y.tree(i).time;

        jT_x = max([state_x.tree(state_x.tree(i).child).time]);
        jT_y = max([state_y.tree(state_y.tree(i).child).time]);

        if i == state_x.root
            assertTrue(testCase, ...
                       mean([iT_x, jT_x]) <= newage_x ...
                           && newage_x <= 2 * iT_x - jT_x);
            assertEqual(testCase, ...
                        logq_x, log(iT_x - jT_x) - log(newage_x - jT_x), ...
                        'AbsTol', 1e-12)

            assertTrue(testCase, ...
                       mean([iT_y, jT_y]) <= newage_y ...
                           && newage_y <= 2 * iT_y - jT_y);
            assertEqual(testCase, ...
                        logq_y, log(iT_y - jT_y) - log(newage_y - jT_y), ...
                        'AbsTol', 1e-12);
        else
            kT_x = state_x.tree(state_x.tree(i).parent).time;
            kT_y = state_y.tree(state_y.tree(i).parent).time;

            assertTrue(testCase, jT_x <= iT_x && iT_x <= kT_x);
            assertTrue(testCase, jT_y <= iT_y && iT_y <= kT_y);

            assertTrue(testCase, jT_x <= newage_x && newage_x <= kT_x);
            assertTrue(testCase, jT_y <= newage_y && newage_y <= kT_y);

            assertEqual(testCase, logq_x, 0);
            assertEqual(testCase, logq_y, 0);
        end
    end

end

function remainCoupledTests(testCase)
    [state, ~] = dummyStates();

    visited = zeros(1, state.NS - 1);
    while all(visited < 1e2)
        [i, newage_x, newage_y, logq_x, logq_y] ...
            = SchooseCoupledMaximal(state, state);
        visited(state.nodes == i) = visited(state.nodes == i) + 1;

        assertEqual(testCase, newage_x, newage_y);
        assertEqual(testCase, logq_x, logq_y);
    end
end

function marginalTests(testCase)
    [state_x, state_y] = dummyStates();

    n = 1e3;
    [xCoupled, yCoupled, xMarginal, yMarginal] = deal(zeros(n, 2));

    for samp = 1:n
        [iCoupled, newage_xCoupled, newage_yCoupled, ~, ~] ...
            = SchooseCoupledMaximal(state_x, state_y);
        xCoupled(samp, :) = [iCoupled, newage_xCoupled];
        yCoupled(samp, :) = [iCoupled, newage_yCoupled];

        [iMarginal_x, newage_xMarginal, ~] = Schoose(state_x);
        xMarginal(samp, :) = [iMarginal_x, newage_xMarginal];

        [iMarginal_y, newage_yMarginal, ~] = Schoose(state_y);
        yMarginal(samp, :) = [iMarginal_y, newage_yMarginal];
    end

    for j = 1:length(state_x.nodes)
        i = state_x.nodes(j);
        subplot(state_x.NS - 1, 2, 2 * j - 1);
        [n_xC, e_xC] = histcounts(xCoupled(xCoupled(:, 1) == i, 2), 20, ...
                                  'Normalization', 'cdf');
        [n_xM, e_xM] = histcounts(xMarginal(xMarginal(:, 1) == i, 2), 20, ...
                                  'Normalization', 'cdf');
        plot([e_xC; e_xM]', [zeros(2, 1), [n_xC; n_xM]]', ':', 'LineWidth', 2);

        subplot(state_x.NS - 1, 2, 2 * j);
        [n_yC, e_yC] = histcounts(yCoupled(yCoupled(:, 1) == i, 2), 20, ...
                                  'Normalization', 'cdf');
        [n_yM, e_yM] = histcounts(yMarginal(yMarginal(:, 1) == i, 2), 20, ...
                                  'Normalization', 'cdf');
        plot([e_yC; e_yM]', [zeros(2, 1), [n_yC; n_yM]]', ':', 'LineWidth', 2);
    end
    v = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end

function coupledTests(testCase)
    [state_x, state_y] = dummyStates();

    assertEqual(testCase, state_x.root, state_y.root);
    assertEqual(testCase, state_x.nodes, state_y.nodes);
    assertEqual(testCase, state_x.NS - 1, length(state_x.nodes));

    [visited, coupled] = deal(zeros(1, state_x.NS - 1));

    while all(visited < 5e3)
        [i, newage_x, newage_y, ~, ~] = SchooseCoupledMaximal(state_x, state_y);
        ind = find(state_x.nodes == i);
        visited(ind) = visited(ind) + 1;
        if newage_x == newage_y
            coupled(ind) = coupled(ind) + 1;
        end
    end
    overlapObs = coupled ./ visited;
    overlapExp = getOverlap(state_x, state_y);

    fprintf('obs     exp     diff\n');
    fmt = '%-.4f  %-.4f  %+-.4f\n';
    fprintf(fmt, [overlapObs; overlapExp; overlapObs - overlapExp]);
    v = input(['Do these proportions of coupled samples match? ', ...
               'Reply 1 for yes... ']);
    assertEqual(testCase, v, 1);
end

function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    global BORROWING MCMCCAT
    testCase.TestData.BORROWING = BORROWING;
    testCase.TestData.MCMCCAT = MCMCCAT;
end

function teardownOnce(testCase)
    global BORROWING MCMCCAT
    BORROWING = testCase.TestData.BORROWING;
    MCMCCAT = testCase.TestData.MCMCCAT;
end

function [state_x, state_y] = dummyStates()
    L = 7;
    t = 0.5;
    state_x = dummyState(ExpTree(L, t));
    state_y = housekeeping(state_x, dummyState(ExpTree(L, t)));
end

function state = dummyState(s)
    % Same as makestate but without data and associated calculations
    % Assume catastrophes are already on tree
    state = tree2state(s);
    state.claderoot = [];
    state.cat = cellfun('length', {state.tree.catloc});
    state.cat = state.cat(:);
    state.ncat = sum(state.cat);
    state.length = TreeLength(state.tree, state.root);
    % Added because housekeeping now uses MarkRCurs when BORROWING = 0
    state.kappa = rand;
end


function overlap = getOverlap(state_x, state_y)
    overlap = nan(1, length(state_x.nodes));
    for ind = 1:length(state_x.nodes)
        i = state_x.nodes(ind);
        [~, kT_x, jT_x, a_x, b_x] ...
            = SchooseCoupledMaximal.nodeTimesAndRanges(i, state_x);
        [~, kT_y, jT_y, a_y, b_y] ...
            = SchooseCoupledMaximal.nodeTimesAndRanges(i, state_y);

        if i == state_x.root
            overlap(ind) = max(0, min(b_x, b_y) - max(a_x, a_y)) ...
                           / max(b_x - a_x, b_y - a_y);
        else
            overlap(ind) = max(0, min(kT_x, kT_y) - max(jT_x, jT_y)) ...
                           / max(kT_x - jT_x, kT_y - jT_y);
        end
    end
end

function tests = sampleCouplingTest
    % Compare to sampleMarginal
    tests = functiontests(localfunctions);
end

% Checking two trees stay together
function cladesNoMoveNarrowStayCoupledTest(testCase)
    moveNarrowStayCoupled(testCase, 'No');
end

function cladesYesMoveNarrowStayCoupledTest(testCase)
    moveNarrowStayCoupled(testCase, 'Yes');
end

function moveNarrowStayCoupled(testCase, clades)
    [s, nReps, newageCoupling, newageMarginal] = getTree10Params(clades);
    nRep = 1;
    while nRep <= nReps
        i = sampleSubtree(s);
        [j, k, FAIL] = BchooseCoupled.getNarrowDestination(i, s);
        if ~FAIL
            THETA = BchooseCoupled.sampleCoupling.sampleTheta();
            [newage_x, newage_y, logq_x, logq_y] ...
                = BchooseCoupled.sampleCoupling(i, j, j, k, k, s, s, THETA);
            assertEqual(testCase, newage_x, newage_y);
            assertEqual(testCase, logq_x, logq_y);
            newageCoupling(nRep) = newage_x;
            [newageMarginal(nRep), ~] = BchooseCoupled.sampleMarginal(...
                i, j, k, s, THETA);
            nRep = nRep + 1;
        end
    end
    BchooseCoupled.sampleCoupling.compareDistributions(testCase, ...
        newageCoupling, newageMarginal);
end

function cladesNoMoveWideStayCoupledTest(testCase)
    [s, nReps, newageCoupling, newageMarginal] = getTree10Params('No');
    nRep = 1;
    while nRep <= nReps
        i = sampleSubtree(s);
        r = 1:(length(s) - 1);
        if length(r) > 4
            [j, ~, k, ~, FAIL, ~] = BchooseCoupled.getWideDestination(i, r, s, s);
        else
            FAIL = 1;
        end
        if ~FAIL
            THETA = BchooseCoupled.sampleCoupling.sampleTheta();
            [newage_x, newage_y, logq_x, logq_y] ...
                = BchooseCoupled.sampleCoupling(i, j, j, k, k, s, s, THETA);
            assertEqual(testCase, newage_x, newage_y);
            assertEqual(testCase, logq_x, logq_y);
            nRep = nRep + 1;
            newageCoupling(nRep) = newage_x;
            [newageMarginal(nRep), ~] = BchooseCoupled.sampleMarginal(...
                i, j, k, s, THETA);
            nRep = nRep + 1;
        end
    end
    BchooseCoupled.sampleCoupling.compareDistributions(testCase, ...
        newageCoupling, newageMarginal);
end

function cladesYesMoveWideStayCoupledTest(testCase)
    [s, nReps, newageCoupling, newageMarginal] = getTree10Params('Yes');
    nRep = 1;
    while nRep <= nReps
        i = sampleSubtree(s);
        r = BchooseCoupled.getWideCandidatesClade(i, s);
        if length(r) > 4
            [j, ~, k, ~, FAIL, ~] = BchooseCoupled.getWideDestination(i, r, s, s);
        else
            FAIL = 1;
        end
        if ~FAIL
            THETA = BchooseCoupled.sampleCoupling.sampleTheta();
            [newage_x, newage_y, logq_x, logq_y] ...
                = BchooseCoupled.sampleCoupling(i, j, j, k, k, s, s, THETA);
            assertEqual(testCase, newage_x, newage_y);
            assertEqual(testCase, logq_x, logq_y);
            newageCoupling(nRep) = newage_x;
            [newageMarginal(nRep), ~] = BchooseCoupled.sampleMarginal(...
                i, j, k, s, THETA);
            nRep = nRep + 1;
        end
    end
    BchooseCoupled.sampleCoupling.compareDistributions(testCase, ...
        newageCoupling, newageMarginal);
end

% Small coupling examples
function narrowCouplingTest(testCase)
    % We will try to couple node times for two trees (details in getTrees5Params)
    [s_x, s_y, nRep, nReps, matchCount, matchAttempt] = getTrees5Params();
    [logq_xObs, logq_yObs, logq_xExp, logq_yExp] = deal(nan(nReps, 1));
    while nRep <= nReps
        i = sampleSubtree(s_x);
        [j_x, k_x, FAIL_x] = BchooseCoupled.getNarrowDestination(i, s_x);
        [j_y, k_y, FAIL_y] = BchooseCoupled.getNarrowDestination(i, s_y);
        if ~(FAIL_x || FAIL_y)
            [newage_x, newage_y, logq_xObs(nRep), logq_yObs(nRep)] ...
                = BchooseCoupled.sampleCoupling(...
                    i, j_x, j_y, k_x, k_y, s_x, s_y, []);
            if newage_x == newage_y
                matchCount(i) = matchCount(i) + 1;
            end
            matchAttempt(i) = matchAttempt(i) + 1;
            logq_xExp(nRep) = narrowLogQ(s_x, i, j_x, k_x);
            logq_yExp(nRep) = narrowLogQ(s_y, i, j_y, k_y);
            nRep = nRep + 1;
        end
    end
    assertEqual(testCase, logq_xObs, logq_xExp, 'AbsTol', 1e-12);
    assertEqual(testCase, logq_yObs, logq_yExp, 'AbsTol', 1e-12);
    % Coupled move fails in x when i = 1, 6, 7, 10 and y when i = 6, 7, 8, 10
    assertEqual(testCase, find(matchAttempt), [2:5, 9]);
    % With narrow moves, only when i = 4, 5 do the time ranges overlap for
    % coupling to be possible
    assertEqual(testCase, find(matchCount), 4:5);
    % We have that pa(pa(4, 5)) = pa(9) = 8 in both; so when i = 4, 5
    %   * in x we draw t_9 uniformly on [t_3, t_8] = [0, 7.9609]
    %   * in y we draw t_9 uniformly on [t_1, t_8] = [0, 46.342]
    % With a maximal coupling, we propose the same node time with probability
    %   [t_8(x) - t_3(x)] / [t_8(y) - t_1(y)] = t_8(x) / t_8(y) = 0.17179
    propObsNum = sum(matchCount(4:5));
    propObsDen = sum(matchAttempt(4:5));
    propObs = propObsNum / propObsDen;
    propExp = s_x(8).time / s_y(8).time;
    fprintf('Proportion coupled\nObserved %g = %d / %d\nExpected %g\n', ...
            propObs, propObsNum, propObsDen, propExp);
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end

function logq = narrowLogQ(s, i, j, k)
    logq = log(s(k).time - max(s(i).time, s(j).time)) ...
           - log(s(k).time - max(s(s(s(i).parent).child).time));
end

function wideCouplingCladesNoTest(testCase)
    % We will try to couple node times for two trees in getTrees5Params
    THETA = 0.01;
    [s_x, s_y, nRep, nReps, matchCount, matchAttempt] = getTrees5Params();
    N = length(s_x) - 1;
    assertGreaterThan(testCase, N, 4);
    while nRep <= nReps
        i = sampleSubtree(s_x);
        r = randperm(N);
        [j_x, j_y, k_x, k_y, FAIL_x, FAIL_y] ...
            = BchooseCoupled.getWideDestination(i, r, s_x, s_y);
        if ~(FAIL_x || FAIL_y)
            [newage_x, newage_y, ~, ~] = BchooseCoupled.sampleCoupling(...
                    i, j_x, j_y, k_x, k_y, s_x, s_y, THETA);
            if newage_x == newage_y
                matchCount(i) = matchCount(i) + 1;
            end
            matchAttempt(i) = matchAttempt(i) + 1;
            nRep = nRep + 1;
        end
    end
    % Coupled move fails when i = 6 (not selected), 7 or 8 (parent is 6 in y so
    % cannot regraft back onto Adam-root branch and proposal distributions do
    % not overlap in rest of tree
    assertEqual(testCase, find(matchAttempt), [1:5, 9]);
    % For now, we only look at coupling events when i = 4 or 5
    % TODO: check coupling events for all nodes
    % When i = 4, 5, we only couple if j = 1, 2, 3, 6 (each w.p. 1/6 among
    % non-failing options) and
    %   * j = 1 w.p. t_x(6) / t_y(8)
    %   * j = 2 w.p. t_x(7) / t_y(7)
    %   * j = 3 w.p. t_x(8) / t_y(7)
    %   * j = 6 w.p. exp(-theta(t_y(6) - t_x(6)))
    propObsNum = sum(matchCount(4:5));
    propObsDen = sum(matchAttempt(4:5));
    propObs = propObsNum / propObsDen;
    propExp = sum([s_x(6).time / s_y(8).time, ...
                   s_x(7).time / s_y(7).time, ...
                   s_x(8).time / s_y(7).time, ...
                   exp(-THETA * (s_y(6).time - s_x(6).time))]) / 6;
    fprintf('Proportion coupled\nObserved %g = %d / %d\nExpected %g\n', ...
            propObs, propObsNum, propObsDen, propExp);
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end


% Helper functions
function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function [s, nReps, newageObs, newageExp] = getTree10Params(clades)
    s = BchooseCoupled.state10a(clades);
    nReps = 5e4;
    [newageObs, newageExp] = deal(nan(nReps, 1));
end

function [s_x, s_y, nRep, nReps, matchCount, matchAttempt] = getTrees5Params()
    % Two ExpTree(L = 5, Theta = 0.012118) trees
    % x on the left, y on the right, Adam node is 10 in both cases
    %       ________  2  __
    %      |               7_____
    %    __7   _____  3  __|     |
    %   |  |__8                  |
    %   |     |   __  4  __      6__
    % __6     |__9         9__   |
    %   |        |__  5  __|  8__|
    %   |                     |
    %   |___________  1  _____|
    % Housekeeping has been run on y so 9 is the only internal node with a
    % common subtree are, in addition to the root node 6
    testStates = '+BchooseCoupled/+sampleCoupling/statesCheckCoupling';
    s_x = getfield(load(testStates), 'state_x', 'tree');
    s_y = getfield(load(testStates), 'state_y', 'tree');

    nRep = 1;
    nReps = 5e5;
    [matchCount, matchAttempt] = deal(zeros(size(s_x)));
end

function i = sampleSubtree(s)
    % Get valid i for Bchoose, not Echoose > NARROW, only used for testing
    global ROOT
    N = length(s) - 1;
    i = ceil(N * rand);
    while s(i).type == ROOT
       i = ceil(N * rand);
    end
end

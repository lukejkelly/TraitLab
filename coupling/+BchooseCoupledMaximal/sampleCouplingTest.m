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
    [s, nReps, newageCoupling, newageMarginal] = getParams(clades);
    nRep = 1;
    while nRep <= nReps
        i = BchooseCoupledMaximal.sampleSubtree(s);
        [j, k, FAIL] = BchooseCoupledMaximal.getNarrowDestination(i, s);
        if ~FAIL
            THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
            [newage_x, newage_y, logq_x, logq_y] ...
                = BchooseCoupledMaximal.sampleCoupling(i, j, j, k, k, s, s, ...
                                                       THETA);
            assertEqual(testCase, newage_x, newage_y);
            assertEqual(testCase, logq_x, logq_y);
            newageCoupling(nRep) = newage_x;
            [newageMarginal(nRep), ~] = BchooseCoupledMaximal.sampleMarginal(...
                i, j, k, s, THETA);
            nRep = nRep + 1;
        end
    end
    BchooseCoupledMaximal.sampleCoupling.compareDistributions(testCase, ...
        newageCoupling, newageMarginal);
end

function cladesNoMoveWideStayCoupledTest(testCase)
    [s, nReps, newageCoupling, newageMarginal] = getParams('No');
    nRep = 1;
    while nRep <= nReps
        i = BchooseCoupledMaximal.sampleSubtree(s);
        r = 1:(length(s) - 1);
        [j, k, FAIL] = BchooseCoupledMaximal.getWideDestination(i, r, ...
                                                                length(r), s);
        if ~FAIL
            THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
            [newage_x, newage_y, logq_x, logq_y] ...
                = BchooseCoupledMaximal.sampleCoupling(i, j, j, k, k, s, s, ...
                                                       THETA);
            assertEqual(testCase, newage_x, newage_y);
            assertEqual(testCase, logq_x, logq_y);
            nRep = nRep + 1;
            newageCoupling(nRep) = newage_x;
            [newageMarginal(nRep), ~] = BchooseCoupledMaximal.sampleMarginal(...
                i, j, k, s, THETA);
            nRep = nRep + 1;
        end
    end
    BchooseCoupledMaximal.sampleCoupling.compareDistributions(testCase, ...
        newageCoupling, newageMarginal);
end

function cladesYesMoveWideStayCoupledTest(testCase)
    [s, nReps, newageCoupling, newageMarginal] = getParams('Yes');
    nRep = 1;
    while nRep <= nReps
        i = BchooseCoupledMaximal.sampleSubtree(s);
        r = BchooseCoupledMaximal.getWideCandidatesClade(i, s);
        [j, k, FAIL] = BchooseCoupledMaximal.getWideDestination(i, r, ...
                                                                length(r), s);
        if ~FAIL
            THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
            [newage_x, newage_y, logq_x, logq_y] ...
                = BchooseCoupledMaximal.sampleCoupling(i, j, j, k, k, s, s, ...
                                                       THETA);
            assertEqual(testCase, newage_x, newage_y);
            assertEqual(testCase, logq_x, logq_y);
            newageCoupling(nRep) = newage_x;
            [newageMarginal(nRep), ~] = BchooseCoupledMaximal.sampleMarginal(...
                i, j, k, s, THETA);
            nRep = nRep + 1;
        end
    end
    BchooseCoupledMaximal.sampleCoupling.compareDistributions(testCase, ...
        newageCoupling, newageMarginal);
end

% Small coupling examples
function narrowCouplingTest(testCase)
    % We will try to couple node times for two random trees
    [s_x, s_y] = getTrees();

    [matchCount_x, matchCount_y, matchAttempt_x, matchAttempt_y] ...
        = deal(zeros(size(s_x)));
    nRep = 1;
    nReps = 1e3;
    while nRep <= nReps
        i = BchooseCoupledMaximal.sampleSubtree(s_x);
        [j_x, k_x, FAIL_x] = BchooseCoupledMaximal.getNarrowDestination(i, s_x);
        [j_y, k_y, FAIL_y] = BchooseCoupledMaximal.getNarrowDestination(i, s_y);
        if ~(FAIL_x || FAIL_y)
            [newage_x, newage_y, logq_xObs, logq_yObs] ...
                = BchooseCoupledMaximal.sampleCoupling(...
                    i, j_x, j_y, k_x, k_y, s_x, s_y, []);
            if newage_x == newage_y
                matchCount_x(s_x(i).parent) = matchCount_x(s_x(i).parent) + 1;
                matchCount_y(s_y(i).parent) = matchCount_y(s_y(i).parent) + 1;
            end
            matchAttempt_x(s_x(i).parent) = matchAttempt_x(s_x(i).parent) + 1;
            matchAttempt_y(s_y(i).parent) = matchAttempt_y(s_y(i).parent) + 1;

            logq_xExp = log(s_x(k_x).time - max(s_x(i).time, s_x(j_x).time)) ...
                        - log(s_x(k_x).time ...
                              - max(s_x(s_x(s_x(i).parent).child).time));
            logq_yExp = log(s_y(k_y).time - max(s_y(i).time, s_y(j_y).time)) ...
                        - log(s_y(k_y).time ...
                              - max(s_y(s_y(s_y(i).parent).child).time));
            assertEqual(testCase, logq_xObs, logq_xExp, 'AbsTol', 1e-12);
            assertEqual(testCase, logq_yObs, logq_yExp, 'AbsTol', 1e-12);
            nRep = nRep + 1;
        end
    end
end


% % Checking marginal and coupling distributions
% function cladesNoMoveNarrowDistributionTest(testCase)
%     moveNarrowDistribution(testCase, 'No');
% end
%
% function cladesYesMoveNarrowDistributionTest(testCase)
%     moveNarrowDistribution(testCase, 'Yes');
% end
%
% function moveNarrowDistribution(testCase, clades)
%     [s_x, nReps, newageObs, newageExp] = getParams(clades);
%     %
%
%
%     nRep = 1;
%     while nRep <= nReps
%         i = BchooseCoupledMaximal.sampleSubtree(s);
%         [j, k, FAIL] = BchooseCoupledMaximal.getNarrowDestination(i, s);
%         if ~FAIL
%             THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
%             [newage_x, newage_y, logq_x, logq_y] ...
%                 = BchooseCoupledMaximal.sampleCoupling(i, j, j, k, k, s, s, ...
%                                                        THETA);
%             assertEqual(testCase, newage_x, newage_y);
%             assertEqual(testCase, logq_x, logq_y);
%             nRep = nRep + 1;
%         end
%     end
% end
%


% Help functions
function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function [s, nReps, newageObs, newageExp] = getParams(clades)
    s = BchooseCoupledMaximal.state10(clades);
    nReps = 1e3;
    [newageObs, newageExp] = deal(nan(nReps, 1));
end

function [s_x, s_y] = getTrees()
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
    testStates = '+BchooseCoupledMaximal/+sampleCoupling/statesCheckCoupling';
    s_x = getfield(load(testStates), 'state_x', 'tree');
    s_y = getfield(load(testStates), 'state_y', 'tree');
end

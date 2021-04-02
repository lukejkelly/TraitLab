function tests = sampleMarginalTest
    % Compare to core/Bchoose.m
    tests = functiontests(localfunctions);
end

function cladesNoMoveNarrowTest(testCase)
    narrowCheck(testCase, 'No');
end

function cladesYesMoveNarrowTest(testCase)
    narrowCheck(testCase, 'Yes');
end

function cladesNoMoveWideBelowRootTest(testCase)
    wideBelowRootCheck(testCase, 'No');
end

function cladesYesMoveWideBelowRootTest(testCase)
    wideBelowRootCheck(testCase, 'Yes');
end

function cladesNoMoveWideAboveRootTest(testCase)
    wideAboveRootCheck(testCase, 'No');
end

function cladesYesMoveWideAboveRootTest(testCase)
    wideAboveRootCheck(testCase, 'Yes');
end

function narrowCheck(testCase, clades)
    % For now we look at overall sampling behaviour rather than comparing
    % distributions at each branch
    [s, state, mt, prior, nReps, newageObs, newageExp] = getParams(clades, ...
                                                                   'NARROW');
    r = 1;
    while r <= nReps
        THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
        [i, j, k, newageExp_r, logqExp] = Bchoose(state, mt, THETA, prior);
        if isempty(newageExp_r)
            assertEqual(testCase, logqExp, -Inf);
        else
            newageExp(r) = newageExp_r;
            [newageObs(r), logqObs] = BchooseCoupledMaximal.sampleMarginal(...
                    i, j, k, s, THETA);
            assertEqual(testCase, logqObs, logqExp, 'AbsTol', 1e-12);
            r = r + 1;
        end
    end
    BchooseCoupledMaximal.sampleCoupling.compareDistributions(testCase, ...
            newageObs, newageExp);
end

function wideBelowRootCheck(testCase, clades)
    % For now we look at overall sampling behaviour below root
    global ROOT
    [s, state, mt, prior, nReps, newageObs, newageExp] = getParams(clades, ...
                                                                   'WIDE');
    r = 1;
    while r <= nReps
        THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
        [i, j, k, newageExp_r, logqExp] = Bchoose(state, mt, THETA, prior);
        if isempty(newageExp_r)
            assertEqual(testCase, logqExp, -Inf);
        elseif s(j).type ~= ROOT
            newageExp(r) = newageExp_r;
            [newageObs(r), logqObs] = BchooseCoupledMaximal.sampleMarginal(...
                    i, j, k, s, THETA);
            assertEqual(testCase, logqObs, logqExp, 'AbsTol', 1e-12);
            r = r + 1;
        end
    end
    BchooseCoupledMaximal.sampleCoupling.compareDistributions(testCase, ...
            newageObs, newageExp);
end

function wideAboveRootCheck(testCase, clades)
    % Only look at Exp(THETA) samples above root
    global ROOT
    [s, state, mt, prior, nReps, newageObs, newageExp] = getParams(clades, ...
                                                                   'WIDE');
    r = 1;
    while r <= nReps
        THETA = BchooseCoupledMaximal.sampleCoupling.sampleTheta();
        [i, j, k, newageExp_r, logqExp] = Bchoose(state, mt, THETA, prior);
        if isempty(newageExp_r)
            assertEqual(testCase, logqExp, -Inf);
        elseif s(j).type == ROOT
            newageExp(r) = newageExp_r;
            [newageObs(r), logqObs] = BchooseCoupledMaximal.sampleMarginal(...
                    i, j, k, s, THETA);
            assertEqual(testCase, logqObs -logqExp, ...
                        THETA * (newageObs(r) - newageExp(r)), 'AbsTol', 1e-12);
            r = r + 1;
        end
    end
    BchooseCoupledMaximal.sampleCoupling.compareDistributions(testCase, ...
            newageObs, newageExp);
end


function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function [s, state, mt, prior, nReps, newageObs, newageExp] = getParams(...
        clades, moveType)
    global NARROW WIDE
    s = BchooseCoupledMaximal.state10a(clades);
    state.tree = s;
    state.NS = length(state.tree) / 2;
    switch moveType
    case 'NARROW'
        mt = NARROW;
    case 'WIDE'
        mt = WIDE;
    otherwise
        error('moveType must be ''NARROW'' or ''WIDE''');
    end
    prior.isclade = strcmp(clades, 'yes');
    nReps = 5e4;
    [newageObs, newageExp] = deal(nan(nReps, 1));
end

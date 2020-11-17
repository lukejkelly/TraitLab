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

function cladesNoMoveWideTest(testCase)
    wideCheck(testCase, 'No');
end

function cladesYesMoveWideTest(testCase)
    wideCheck(testCase, 'Yes');
end

function narrowCheck(testCase, clades)
    % For now we look at overall sampling behaviour rather than comparing
    % distributions at each branch
    [s, state, mt, THETA, prior, nReps, newageObs, newageExp] ...
            = getParams(clades, 'NARROW');
    r = 1;
    while r <= nReps
        [i, j, k, newageExp_r, logqExp] = Bchoose(state, mt, THETA, prior);
        if isempty(newageExp_r)
            assertEqual(testCase, logqExp, -Inf);
        else
            newageExp(r) = newageExp_r;
            [newageObs(r), logqObs] = BchooseCoupledMaximal.sampleMarginal(...
                    i, j, k, s, THETA);
            assertEqual(testCase, logqObs, logqExp, 'AbsTol', 1e-12);
        end
        r = r + 1;
    end
    compareDistributions(testCase, newageObs, newageExp);
end

function wideCheck(testCase, clades)
    % For now we look at overall sampling behaviour rather than comparing
    % distributions at each branch
    warning('Still some checking to do depending on the value of THETA');
    global ROOT
    [s, state, mt, THETA, prior, nReps, newageObs, newageExp] ...
            = getParams(clades, 'WIDE');
    r = 1;
    while r <= nReps
        [i, j, k, newageExp_r, logqExp] = Bchoose(state, mt, THETA, prior);
        if isempty(newageExp_r)
            assertEqual(testCase, logqExp, -Inf);
        else
            newageExp(r) = newageExp_r;
            [newageObs(r), logqObs] = BchooseCoupledMaximal.sampleMarginal(...
                    i, j, k, s, THETA);
            if s(j).type == ROOT
                jT = s(j).time;
                deltaObs = newageObs(r) - jT;
                deltaExp = newageExp(r) - jT;
                assertEqual(testCase, logqObs -logqExp, ...
                            THETA * (deltaObs - deltaExp), 'AbsTol', 1e-12);
            else
                assertEqual(testCase, logqObs, logqExp, 'AbsTol', 1e-12);
            end
        end
        r = r + 1;
    end
    compareDistributions(testCase, newageObs, newageExp);
end


function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function [s, state, mt, THETA, prior, nReps, newageObs, newageExp] ...
        = getParams(clades, moveType)
    global NARROW WIDE
    s = BchooseCoupledMaximal.state10(clades);
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
    THETA = 5e-4 + rand * 1.5e-3;
    prior.isclade = strcmp(clades, 'yes');
    nReps = 5e3;
    [newageObs, newageExp] = deal(nan(nReps, 1));
end

function compareDistributions(testCase, xObs, xExp)
    [nObs, eObs] = histcounts(xObs, 20, 'Normalization', 'cdf');
    [nExp, eExp] = histcounts(xExp, 20, 'Normalization', 'cdf');
    plot([eObs; eExp]', [zeros(2, 1), [nObs; nExp]]', ':', 'LineWidth', 2);
    legend('Obs', 'Exp');
    v = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end

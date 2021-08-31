function tests = SchooseTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Only checking outputs as distributions checked elsewhere
    global BORROWING MCMCCAT
    for L = 8:12
        for BORROWING = 0:1
            for MCMCCAT = 0:1
                state = dummyState(L);
                rngCurr = rng;
                [iObs, newageObs, logqObs, catObs, locObs] = Schoose(state);
                rng(rngCurr);
                [iExp, newageExp, logq_t] = SchooseTime(state);
                if MCMCCAT
                    [catExp, locExp, logq_c] = SchooseCats(state, iExp, newageExp);
                else
                    catExp = [];
                    locExp = [];
                    logq_c = 0;
                end
                rng('shuffle');
                logqExp = logq_t + logq_c;

                assertEqual(testCase, [iObs, newageObs, logqObs], ...
                            [iExp, newageExp, logqExp]);
                assertEqual(testCase, catObs, catExp);
                assertEqual(testCase, locObs, locExp);
            end
        end
    end
end

function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function state = dummyState(L)
    global MCMCCAT
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    if MCMCCAT
        for i = 1:(5 + poissrnd(5))
            state = AddCat(state);
        end
    end
end

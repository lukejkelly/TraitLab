function tests = LogPriorParmTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    global VARYNU VARYBETA

    state = struct('mu', rand, 'beta', rand);
    nstate = struct('mu', rand, 'beta', rand);

    VARYNU = 0;
    VARYBETA = 0;
    obs1 = LogPriorParm(state, nstate);
    VARYBETA = 1;
    obs2 = LogPriorParm(state, nstate);


    exp1 = -0.999 * (log(nstate.mu) - log(state.mu)) ...
            - 0.001 * (nstate.mu - state.mu);
    exp2 = exp1 + (-0.999 * (log(nstate.beta) - log(state.beta)) ...
             - 0.001 * (nstate.beta - state.beta));

    assertEqual(testCase, obs1, exp1, 'AbsTol', 1e-6);
    assertEqual(testCase, obs2, exp2, 'AbsTol', 1e-6);
end

function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    testCase.TestData.VARYNU = VARYNU;
    testCase.TestData.VARYBETA = VARYBETA;
end

function teardownOnce(testCase)
    global VARYNU VARYBETA
    VARYNU = testCase.TestData.VARYNU;
    VARYBETA = testCase.TestData.VARYBETA;
end

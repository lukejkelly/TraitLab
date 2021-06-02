function tests = scaleCoupledTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    global DONTMOVECATS;
    mcmc = struct('update', struct('del', 0.5, 'deldel', 1.5));
    testCase.TestData.mcmc = mcmc;
    testCase.TestData.n_i = 2e1;
    testCase.TestData.n_j = 1e4;
    testCase.TestData.DONTMOVECATS = DONTMOVECATS;
    DONTMOVECATS = 0;
end

function teardownOnce(testCase)
    global DONTMOVECATS
    DONTMOVECATS = testCase.TestData.DONTMOVECATS;
    close;
end

function allTest(testCase)
    allTestRoutine(testCase, 'beta', @scaleCoupledBeta, ...
                   @scaleCoupledTest.scaleCoupledBeta);
    allTestRoutine(testCase, 'kappa', @scaleCoupledKappa, ...
                   @scaleCoupledTest.scaleCoupledKappa);
    allTestRoutine(testCase, 'lambda', @scaleCoupledLambda, ...
                   @scaleCoupledTest.scaleCoupledLambda);
    allTestRoutine(testCase, 'mu', @scaleCoupledMu, ...
                   @scaleCoupledTest.scaleCoupledMu);
    allTestRoutine(testCase, 'rho', @scaleCoupledRho, ...
                   @scaleCoupledTest.scaleCoupledRho);
    warning('Write test routine for coupling xi');

    title('Coupling proportions');
    fprintf('Figure shows proportion of matching samples\n');
    fprintf('%i trials of %i samples\n', testCase.TestData.n_i, ...
            testCase.TestData.n_j);
    v = input('Are these proportions okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function allTestRoutine(testCase, par, testFun, legacyFun)
    coupledTestRoutine(testCase, par, testFun);
    legacyTestRoutine(testCase, par, testFun, legacyFun);
end

function coupledTestRoutine(testCase, par, testFun)
    n_i = testCase.TestData.n_i;
    n_j = testCase.TestData.n_j;
    mcmc = testCase.TestData.mcmc;

    o = nan(n_i, 1);
    c = nan(n_i, n_j);
    for i = 1:n_i
        state_x = scaleCoupledTest.dummyState();
        state_y = scaleCoupledTest.dummyState();
        o(i) = scaleCoupledTest.getOverlap(state_x.(par), state_y.(par), mcmc);
        for j = 1:n_j
            [nstate_x, nstate_y, ~, ~, ~, ~, logq_x, logq_y] ...
                = testFun(state_x, state_y, mcmc);
            c(i, j) = ismembertol(nstate_x.(par), nstate_y.(par));
        end
        % Check once per data pair of states
        assertEqual(testCase, nstate_x.(par) / state_x.(par), ...
                    exp(-logq_x), 'RelTol', 1e-12);
        assertEqual(testCase, nstate_y.(par) / state_y.(par), ...
                    exp(-logq_y), 'RelTol', 1e-12);
    end
    % Proportion of matching samples
    subplot(3, 2, find(cellfun(@(z) strcmp(par, z) , ...
                               {'beta', 'kappa', 'lambda', 'mu', 'rho'})));
    plot(o, mean(c, 2) - o, 'x');
    subtitle(par);
    xlabel('Observed');
    ylabel('Difference');
    refline(0, 0);
end

function legacyTestRoutine(testCase, par, testFun, legacyFun)
    % compare new and old versions
    n_i = testCase.TestData.n_i;
    mcmc = testCase.TestData.mcmc;
    for i = 1:n_i
        state_x = scaleCoupledTest.dummyState();
        state_y = scaleCoupledTest.dummyState();

        rng_state = rng;
        [nstate_x1, nstate_y1, U_x1, U_y1, OK_x1, OK_y1, logq_x1, logq_y1] ...
            = testFun(state_x, state_y, mcmc);
        rng(rng_state);
        [nstate_x2, nstate_y2, U_x2, U_y2, OK_x2, OK_y2, logq_x2, logq_y2] ...
            = legacyFun(state_x, state_y, mcmc);
        rng('shuffle');

        assertEqual(testCase, nstate_x1.(par), nstate_x2.(par));
        assertEqual(testCase, nstate_y1.(par), nstate_y2.(par));

        assertEqual(testCase, U_x1, U_x2);
        assertEqual(testCase, U_y1, U_y2);

        assertEqual(testCase, OK_x1, OK_x2);
        assertEqual(testCase, OK_y1, OK_y2);

        assertEqual(testCase, logq_x1, logq_x2);
        assertEqual(testCase, logq_y1, logq_y2);

        switch par
        case 'kappa'
            % Legacy code only updates if OK
            if OK_x2
                assertEqual(testCase, nstate_x1.nu, nstate_x2.nu);
            else
                assertEqual(testCase, state_x.nu, nstate_x2.nu);
            end
            if OK_y2
                assertEqual(testCase, nstate_y1.nu, nstate_y2.nu);
            else
                assertEqual(testCase, state_y.nu, nstate_y2.nu);
            end
        case {'lambda', 'mu'}
            assertEqual(testCase, nstate_x1.nu, nstate_x2.nu);
            assertEqual(testCase, nstate_y1.nu, nstate_y2.nu);
        end
    end
end

function tests = checkCouplingTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    % Global variables from pausing runmcmcCoupled and entering an otherwise
    % empty workspace and saving what remained
    load('coupling/+MarkovCoupledTest/globalVariables-20210603.mat');
    testCase.TestData.global_vars = whos('global');
end

function setup(~)
    clear logLkd2 patternCounts;
end

function teardownOnce(testCase)
    cellfun(@clear, {testCase.TestData.global_vars.name});
end

function teardown(~)
    clear logLkd2 patternCounts;
end

function coupledTest(testCase)
    global BORROWING
    load('coupling/+MarkovCoupledTest/coupledVariables-20210603.mat', ...
         'mcmc', 'model', 'state_x', 'state_y');
    assertTrue(testCase, checkCoupling(state_x, state_x));
    assertTrue(testCase, checkCoupling(state_x, state_y));
    assertTrue(testCase, checkCoupling(state_y, state_x));
    assertTrue(testCase, checkCoupling(state_y, state_y));

    % Making small modifications and checking again
    nodes = state_x.nodes(randperm(length(state_x.nodes)));

    % Shift node time
    state_x1 = state_x;
    i1 = nodes(1);
    c1 = max(state_x1.tree(state_x1.tree(i1).child).time);
    p1 = state_x1.tree(state_x1.tree(i1).parent).time;
    state_x1.tree(i1).time = c1 + rand * (p1 - c1);
    state_x1 = updateStateVariables(state_x1);
    assertFalse(testCase, checkCoupling(state_x, state_x1));
    assertFalse(testCase, checkCoupling.extraChecks(state_x, state_x1));
    state_x1h = updateStateVariables(housekeeping(state_x, state_x1));
    assertFalse(testCase, checkCoupling(state_x, state_x1h));
    assertFalse(testCase, checkCoupling.extraChecks(state_x, state_x1h));

    % Switch siblings
    state_x2 = state_x;
    i2 = nodes(2);
    c2 = state_x2.tree(i2).child;
    [state_x2.tree(c2).sibling] = deal(2, 1);
    state_x2.tree(i2).child = fliplr(state_x2.tree(i2).child);
    state_x2 = updateStateVariables(state_x2);
    % This would never happen as housekeeping would occur first
    assertFalse(testCase, checkCoupling(state_x, state_x2));
    assertFalse(testCase, checkCoupling.extraChecks(state_x, state_x2));
    state_x2h = updateStateVariables(housekeeping(state_x, state_x2));
    assertTrue(testCase, checkCoupling(state_x, state_x2h));
    assertTrue(testCase, checkCoupling.extraChecks(state_x, state_x2h));
    
    % Draw samples and compare, catastrophe outputs depend on model
    mcmc.subsample = 5e1;
    for n = 1:5e1
        BORROWING = 0;
        [state_x, ~] = Markov(mcmc, model, state_x);
        assertTrue(testCase, checkCoupling(state_x, state_x));        
    end
    for n = 1:5e1
        BORROWING = 1;
        [state_y, ~] = Markov(mcmc, model, state_y);
        assertTrue(testCase, checkCoupling(state_y, state_y));        
    end
end

function uncoupledTest(testCase)
    global BORROWING
    load('coupling/+MarkovCoupledTest/uncoupledVariables-20210331.mat', ...
         'mcmc', 'model', 'state_x', 'state_y');
    assertFalse(testCase, checkCoupling(state_x, state_y));
    assertFalse(testCase, checkCoupling(state_y, state_x));
    
    % Draw samples and compare, catastrophe outputs depend on model
    mcmc.subsample = 5e1;
    for n = 1:5e1
        BORROWING = 0;
        [state_xs, pa_xs] = Markov(mcmc, model, state_x);
        if any(pa_xs > 0)
            assertFalse(testCase, ...
                        checkCoupling(state_x, ...
                                      housekeeping(state_x, state_xs)));
        else
            assertTrue(testCase, ...
                       checkCoupling(state_x, housekeeping(state_x, state_xs)));
        end
    end
    for n = 1:5e1
        BORROWING = 1;
        [state_ys, pa_ys] = Markov(mcmc, model, state_y);
        if any(pa_ys > 0)
            assertFalse(testCase, ...
                        checkCoupling(state_y, ...
                                      housekeeping(state_y, state_ys)));        
        else
            assertTrue(testCase, ...
                       checkCoupling(state_y, housekeeping(state_y, state_ys)));
        end
    end

end

function state = updateStateVariables(state)
    % TODO: This doesn't change prior field
    global BORROWING
    if BORROWING
        [state.loglkd, state.fullloglkd] = logLkd2(state);
    else
        state.loglkd = LogLkd(state);
        state.fullloglkd = LogLkd(state, state.lambda);
    end
    state.length = TreeLength(state.tree, state.root);

end

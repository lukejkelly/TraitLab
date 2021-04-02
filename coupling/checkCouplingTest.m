function tests = checkCouplingTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    testCase.TestData.global_vars = whos('global');
end

function teardownOnce(testCase)
    cellfun(@clear, {testCase.TestData.global_vars.name});
end

function coupledTest(testCase)
    load('coupling/MarkovTestCoupled-20210330.mat', 'state_x', 'state_y');
    assertTrue(testCase, checkCoupling(state_x, state_x));
    assertTrue(testCase, checkCoupling(state_x, state_y));
    assertTrue(testCase, checkCoupling(state_y, state_x));
    assertTrue(testCase, checkCoupling(state_y, state_y));

    % Making small modifications and checking again
    nodes = state_x.nodes(randperm(length(state_x.nodes)));
    leaves = state_x.leaves(randperm(length(state_x.leaves)));
    
    % Shift node time    
    state_x1 = state_x;
    i1 = nodes(1);
    c1 = max(state_x1.tree(state_x1.tree(i1).child).time);
    p1 = state_x1.tree(state_x1.tree(i1).parent).time;
    state_x1.tree(i1).time = c1 + rand * (p1 - c1);
    state_x1 = updateStateVariables(state_x1);
    assertFalse(testCase, checkCoupling(state_x, state_x1));
    assertFalse(testCase, checkCouplingExtra(state_x, state_x1));
    state_x1h = updateStateVariables(housekeeping(state_x, state_x1));
    assertFalse(testCase, checkCoupling(state_x, state_x1h));
    assertFalse(testCase, checkCouplingExtra(state_x, state_x1h));
   
    % Switch siblings
    state_x2 = state_x;
    i2 = nodes(2);
    c2 = state_x2.tree(i2).child;
    [state_x2.tree(c2).sibling] = deal(2, 1);
    state_x2.tree(i2).child = fliplr(state_x2.tree(i2).child);
    state_x2 = updateStateVariables(state_x2);
    % This would never happen as housekeeping would occur first    
    assertTrue(testCase, checkCoupling(state_x, state_x2));
    assertFalse(testCase, checkCouplingExtra(state_x, state_x2));
    state_x2h = updateStateVariables(housekeeping(state_x, state_x2));
    assertTrue(testCase, checkCoupling(state_x, state_x2h));
    assertTrue(testCase, checkCouplingExtra(state_x, state_x2h));
end

function uncoupledTest(testCase)
    load('coupling/MarkovTestUncoupled-20210331.mat', 'state_x', 'state_y');
    assertFalse(testCase, checkCoupling(state_x, state_y));
    assertFalse(testCase, checkCoupling(state_y, state_x));
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
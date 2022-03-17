function tests = checkCouplingTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    load('coupling/+MarkovCoupledTest/coupledVariables-20210803.mat', ...
         'mcmc', 'model', 'state_x', 'state_y');
    state_x.logprior = LogPrior(model.prior, state_x);
    state_y.logprior = LogPrior(model.prior, state_y);
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
    if i1 == state_x1.root
        t1 = state_x1.tree(i1).time;
        state_x1.tree(i1).time = (c1 + t1) / 2 + rand * (t1 - c1) * 3 / 2;
    else
        p1 = state_x1.tree(state_x1.tree(i1).parent).time;
        state_x1.tree(i1).time = c1 + rand * (p1 - c1);
    end
    state_x1 = updateStateVariables(state_x1, model.prior);
    assertFalse(testCase, checkCoupling(state_x, state_x1));
    assertFalse(testCase, checkCoupling.extraChecks(state_x, state_x1));
    state_x1h = updateStateVariables(housekeeping(state_x, state_x1), model.prior);
    assertFalse(testCase, checkCoupling(state_x, state_x1h));
    assertFalse(testCase, checkCoupling.extraChecks(state_x, state_x1h));

    % Switch siblings
    state_x2 = state_x;
    i2 = nodes(2);
    c2 = state_x2.tree(i2).child;
    [state_x2.tree(c2).sibling] = deal(2, 1);
    state_x2.tree(i2).child = fliplr(state_x2.tree(i2).child);
    state_x2 = updateStateVariables(state_x2, model.prior);
    % This would never happen as housekeeping would occur first
    assertFalse(testCase, checkCoupling(state_x, state_x2));
    assertFalse(testCase, checkCoupling.extraChecks(state_x, state_x2));
    state_x2h = updateStateVariables(housekeeping(state_x, state_x2), model.prior);
    assertTrue(testCase, checkCoupling(state_x, state_x2h));
    assertTrue(testCase, checkCoupling.extraChecks(state_x, state_x2h));

    % Change one xi
    state_x3 = state_x;
    i3 = state_x3.leaves(ceil(rand * length(state_x3.leaves)));
    state_x3.tree(i3).xi = state_x3.tree(i3).xi * 0.99;
    state_x3 = updateStateVariables(state_x3, model.prior);
    assertFalse(testCase, checkCoupling.compareTrees(state_x, state_x3));
    assertFalse(testCase, checkCoupling(state_x, state_x3));

    % Draw samples and compare, catastrophe outputs depend on model
    mcmc.subsample = 10;
    for n = 1:5e1
        [state_x, ~] = Markov(mcmc, model, state_x);
        assertTrue(testCase, checkCoupling(state_x, state_x));
        [state_y, ~] = Markov(mcmc, model, state_y);
        assertTrue(testCase, checkCoupling(state_y, state_y));
    end
end

%% Assertions in this test don't always hold
% function uncoupledTest(testCase)
%     load('coupling/+MarkovCoupledTest/uncoupledVariables-20210803.mat', ...
%          'mcmc', 'model', 'state_x', 'state_y');
%     state_x.logprior = LogPrior(model.prior, state_x);
%     state_y.logprior = LogPrior(model.prior, state_y);
%     assertFalse(testCase, checkCoupling(state_x, state_y));
%     assertFalse(testCase, checkCoupling(state_y, state_x));
%
%     mcmc.subsample = 10;
%     for n = 1:5e1
%         [state_xs, pa_xs] = Markov(mcmc, model, state_x);
%         state_xh = housekeeping(state_x, state_xs);
%         if any(pa_xs > 0) && any(find(pa_xs > 0) ~= 15)
%             % Move 15 may propose 0 catastrophes where there was 0 already
%             assertFalse(testCase, checkCoupling(state_x, state_xh));
%         else
%             assertTrue(testCase, checkCoupling(state_x, state_xh));
%         end
%     end
%     [state_ys, pa_ys] = Markov(mcmc, model, state_y);
%     state_yh = housekeeping(state_y, state_ys);
%     if any(pa_ys > 0) && any(find(pa_ys) ~= 15)
%         assertFalse(testCase, checkCoupling(state_y, state_yh));
%     else
%         assertTrue(testCase, checkCoupling(state_y, state_yh));
%     end
% end

function catTest(testCase)
    load('coupling/+MarkovCoupledTest/uncoupledVariables-20210803.mat', ...
         'state_x', 'model');
    state_x.logprior = LogPrior(model.prior, state_x);
    while all(state_x.cat < 3)
        state_x = AddCat(state_x);
    end
    cInds = find(state_x.cat);

    assertTrue(testCase, checkCoupling(state_x, state_x));

    % Changing catastrophe count
    state_y1 = state_x;
    i = cInds(ceil(rand * length(cInds)));
    state_y1.ncat = state_y1.ncat + 1;
    state_y1.cat(i) = state_y1.cat(i) + 1;
    state_y1 = updateStateVariables(state_y1, model.prior);

    assertFalse(testCase, checkCoupling.compareStates(state_x, state_y1));
    assertFalse(testCase, checkCoupling(state_x, state_y1));

    state_y2 = state_x;
    i = cInds(ceil(rand * length(cInds)));
    state_y2.ncat = state_y2.ncat - 1;
    state_y2.cat(i) = state_y2.cat(i) - 1;
    state_y2 = updateStateVariables(state_y2, model.prior);

    assertFalse(testCase, checkCoupling.compareStates(state_x, state_y2));
    assertFalse(testCase, checkCoupling(state_x, state_y2));

    % Shift cat location
    state_y3 = state_x;
    i = cInds(ceil(rand * length(cInds)));
    state_y3.tree(i).catloc(1) = rand;
    state_y3 = updateStateVariables(state_y3, model.prior);

    assertFalse(testCase, checkCoupling.compareStates(state_x, state_y3));
    assertFalse(testCase, checkCoupling.compareTrees(state_x, state_y3));
    assertFalse(testCase, checkCoupling(state_x, state_y3));

    % Swap cat locations on branch
    state_y4 = state_x;
    i = find(state_y4.cat == 3, 1);
    state_y4.tree(i).catloc = fliplr(state_y4.tree(i).catloc);
    state_y4 = updateStateVariables(state_y4, model.prior);

    assertFalse(testCase, checkCoupling.compareStates(state_x, state_y4));
    assertFalse(testCase, checkCoupling.compareTrees(state_x, state_y4));
    assertFalse(testCase, checkCoupling(state_x, state_y4));
end


function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    load('coupling/+MarkovCoupledTest/globalVariables-20210803.mat');
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

function state = updateStateVariables(state, prior)
    global BORROWING
    if BORROWING
        [state.loglkd, state.fullloglkd] = logLkd2(state);
    else
        state.loglkd = LogLkd(state);
        state.fullloglkd = LogLkd(state, state.lambda);
    end
    state.logprior = LogPrior(prior, state);
    state.length = TreeLength(state.tree, state.root);
end

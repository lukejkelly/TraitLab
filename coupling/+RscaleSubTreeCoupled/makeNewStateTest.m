function tests = makeNewStateTest
    tests = functiontests(localfunctions);
end

function rootPassTest(testCase)
    % Identical to Rscale
    global LEAF;
    state = testCase.TestData.state;
    i = state.root;
    iT = state.tree(state.root).time;
    nu = progeny(state.tree, i, LEAF);

    [lf, nd, t0] = RscaleSubTreeCoupled.getNodePars(state, nu);
    assertTrue(testCase, all(sort(lf) == sort(state.leaves)));
    assertTrue(testCase, all(sort(nd) == sort(state.nodes)));
    assertEqual(testCase, t0, min([state.tree(nu(1, :)).time]));

    rho = 1.1;
    newage = t0 + rho * (iT - t0);

    % Manual checks
    [nstateObs, UObs, logqObs, OKObs] ...
        = RscaleSubTreeCoupled.makeNewState(state, i, lf, nd, t0, newage);
    assertEqual(testCase, nstateObs.tree(i).time, newage, 'AbsTol', 1e-10);
    for j = state.nodes
        assertEqual(testCase, (nstateObs.tree(j).time - t0) / (newage - t0), ...
                    (state.tree(j).time - t0) / (iT - t0), 'AbsTol', 1e-10);
    end
    assertNotEqual(testCase, nstateObs.length, state.length);
    assertTrue(testCase, all(ismember(UObs, state.nodes)));
    assertEqual(testCase, logqObs, (length(state.nodes) - 2) * log(rho), ...
                'AbsTol', 1e-10);
    assertEqual(testCase, OKObs, 1);

    % Comparison with Rscale
    [nstateExp, UExp, ~, OKExp, logqExp] = Rscale(state, rho);
    for j = state.nodes
        assertEqual(testCase, nstateObs.tree(j).time, ...
                    nstateExp.tree(j).time, 'AbsTol', 1e-10);
    end
    assertNotEqual(testCase, nstateObs.length, state.length);
    assertTrue(testCase, all(unique(UObs) == sort(UExp)));
    assertEqual(testCase, OKObs, OKExp);
    assertEqual(testCase, logqObs, logqExp);
end

function rootFailTest(testCase)
    % Identical to Rscale
    global LEAF;
    state = testCase.TestData.state;
    i = state.root;
    iT = state.tree(state.root).time;
    nu = progeny(state.tree, i, LEAF);

    [lf, nd, t0] = RscaleSubTreeCoupled.getNodePars(state, nu);
    assertTrue(testCase, all(sort(lf) == sort(state.leaves)));
    assertTrue(testCase, all(sort(nd) == sort(state.nodes)));
    assertEqual(testCase, t0, min([state.tree(nu(1, :)).time]));

    rho = 0.51;
    newage = t0 + rho * (iT - t0);

    % Manual checks
    [nstateObs, UObs, logqObs, OKObs] ...
        = RscaleSubTreeCoupled.makeNewState(state, i, lf, nd, t0, newage);
    assertEqual(testCase, nstateObs.tree(i).time, newage, 'AbsTol', 1e-10);
    for j = state.nodes
        assertEqual(testCase, (nstateObs.tree(j).time - t0) / (newage - t0), ...
                    (state.tree(j).time - t0) / (iT - t0), 'AbsTol', 1e-10);
    end
    % Values returned when not OK
    assertEqual(testCase, nstateObs.length, state.length);
    assertEmpty(testCase, UObs);
    assertEqual(testCase, logqObs, 0);
    assertEqual(testCase, OKObs, 0);

    % Comparison with Rscale
    [nstateExp, UExp, ~, OKExp, logqExp] = Rscale(state, rho);
    for j = state.nodes
        assertEqual(testCase, nstateObs.tree(j).time, ...
                    nstateExp.tree(j).time, 'AbsTol', 1e-10);
    end
    % Rscale computes quantities anyway
    assertEqual(testCase, nstateExp.length, state.length);
    assertEqual(testCase, sort(UExp), sort(state.nodes));
    assertEqual(testCase, OKExp, 0);
    assertEqual(testCase, logqExp, (length(nd) - 2) * log(rho));
end

function nodePassTest(testCase)
    % Parent of subtree with leaves '2', '7' and '8'
    global LEAF;
    state = testCase.TestData.state;
    i = 6;
    iT = state.tree(i).time;
    nu = progeny(state.tree, i, LEAF);

    [lf, nd, t0] = RscaleSubTreeCoupled.getNodePars(state, nu);
    assertTrue(testCase, all(sort(lf) == intersect(nu(1, :), state.leaves)));
    assertTrue(testCase, all(sort(nd) == intersect(nu(1, :), state.nodes)));
    assertEqual(testCase, t0, min([state.tree(nu(1, :)).time]));

    rho = 1.2;
    newage = t0 + rho * (iT - t0);

    [nstate, U, logq, OK] = RscaleSubTreeCoupled.makeNewState(state, i, lf, ...
                                                              nd, t0, newage);
    assertEqual(testCase, nstate.tree(i).time, newage, 'AbsTol', 1e-10);
    for j = state.nodes
        if ismember(j, nu(1, :))
            assertEqual(testCase, (nstate.tree(j).time - t0) / (newage - t0), ...
                        (state.tree(j).time - t0) / (iT - t0), 'AbsTol', 1e-10);
        else
            assertEqual(testCase, nstate.tree(j).time, state.tree(j).time);
        end
    end
    assertNotEqual(testCase, nstate.length, state.length);
    assertTrue(testCase, all(sort(U) == sort(above([6, 7], nstate.tree, ...
                                                   nstate.root))));
    % Only scaling 2 nodes so terms cancel
    assertEqual(testCase, logq, 0, 'AbsTol', 1e-10);
    assertEqual(testCase, OK, 1);
end

function nodeFailTest(testCase)
    % Fail when subtree on leaves '6', '9' and '10' root older than parent
    global LEAF;
    state = testCase.TestData.state;
    i = 13;
    iT = state.tree(i).time;
    nu = progeny(state.tree, i, LEAF);

    [lf, nd, t0] = RscaleSubTreeCoupled.getNodePars(state, nu);
    assertTrue(testCase, all(sort(lf) == intersect(nu(1, :), state.leaves)));
    assertTrue(testCase, all(sort(nd) == intersect(nu(1, :), state.nodes)));
    assertEqual(testCase, t0, min([state.tree(nu(1, :)).time]));

    rho = 1.5;
    newage = t0 + rho * (iT - t0);

    [nstate, U, logq, OK] = RscaleSubTreeCoupled.makeNewState(state, i, lf, ...
                                                              nd, t0, newage);
    assertEqual(testCase, nstate.tree(i).time, newage, 'AbsTol', 1e-10);
    for j = state.nodes
        if ismember(j, nu(1, :))
            assertEqual(testCase, (nstate.tree(j).time - t0) / (newage - t0), ...
                        (state.tree(j).time - t0) / (iT - t0), 'AbsTol', 1e-10);
        else
            assertEqual(testCase, nstate.tree(j).time, state.tree(j).time);
        end
    end
    assertEqual(testCase, nstate.length, state.length);
    assertTrue(testCase, isempty(U));
    assertEqual(testCase, logq, 0);
    assertEqual(testCase, OK, 0);
end

% Helper functions
function setupOnce(testCase)
    RscaleSubTreeCoupled.unitTests.setupOnce(testCase);
    % 10 leaf tree
    s = rnextree('(((4:315.995,3:10):499.7965,((8:135.1835,7:135.1835):30.462,2:435.6455):380.1461):28.2761,((((9:23.8604,10:23.8604):194.0395,6:217.8999):27.6933,5:245.5932):321.0927,1:566.6859):277.3816)');
    state = RscaleSubTreeCoupled.unitTests.dummyState(s);
    testCase.TestData.state = state;
end

function tearDownOnce(testCase)
    RscaleSubTreeCoupled.unitTests.tearDownOnce(testCase);
end

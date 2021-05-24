function tests = getNodeParsTest
    tests = functiontests(localfunctions);
end

function overallTest(testCase)
    global LEAF ANST ROOT;
    for L = 5:15
        for i = 1:5
            state = RscaleSubTreeCoupled.unitTests.dummyState(ExpTree(L, 1));
            s = state.tree;
            for j = state.nodes
                nu = progeny(s, j, LEAF);

                [lfObs, ndObs, t0Obs] = RscaleSubTreeCoupled.getNodePars(state, nu);

                lfExp = nu(1, [s(nu(1, :)).type] == LEAF);
                ndExp = nu(1, arrayfun(@(k) any(s(k).type == [ANST, ROOT]), ...
                                       nu(1, :)));
                t0Exp = min([s(lfExp).time]);

                assertTrue(testCase, all(sort(lfObs) == sort(lfExp)));
                assertTrue(testCase, all(sort(ndObs) == sort(ndExp)));
                assertEqual(testCase, t0Obs, t0Exp);

                Vel = nu(1, :);
                assertEqual(testCase, length(Vel) - length(lfObs), length(ndObs));
            end
        end
    end
end

function setupOnce(testCase)
    RscaleSubTreeCoupled.unitTests.setupOnce(testCase);
end

function tearDownOnce(testCase)
    RscaleSubTreeCoupled.unitTests.tearDownOnce(testCase);
end

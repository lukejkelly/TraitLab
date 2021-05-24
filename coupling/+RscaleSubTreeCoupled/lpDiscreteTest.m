function tests = lpDiscreteTest
    tests = functiontests(localfunctions);
end

function overallTest(testCase)
    global LEAF;
    for L = 5:15
        for i = 1:5
            state = RscaleSubTreeCoupled.unitTests.dummyState(ExpTree(L, 1));
            for j = state.nodes
                nu = progeny(state.tree, j, LEAF);
                lp = @(x) RscaleSubTreeCoupled.lpDiscrete(x, nu);

                if j == state.root
                    lpr = arrayfun(lp, state.nodes);
                    assertTrue(testCase, all(lpr < 0));
                    assertEqual(testCase, sum(exp(lpr)), 1, 'RelTol', 1e-6);
                    lpl = arrayfun(lp, state.leaves);
                    assertTrue(testCase, all(isinf(lpl)));
                end
                nu0 = nu(1, nu(2, :) == 0);
                lp0 = arrayfun(lp, nu0);
                assertTrue(testCase, all(isinf(lp0)));

                nu1 = nu(1, nu(2, :) ~= 0);
                lp1 = arrayfun(lp, nu1);
                if length(nu1) == 1
                    assertEqual(testCase, lp1, 0);
                else
                    assertTrue(testCase, all(lp1 < 0));
                end
                assertEqual(testCase, sum(exp(lp1)), 1, 'RelTol', 1e-6);
                for k = 1:size(nu, 2)
                   assertEqual(testCase, lp(nu(1, k)), ...
                               log(nu(2, k)) - log(sum(nu(2, :))), ...
                               'AbsTol', 1e-6);
                end
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

function tests = getAgeParsTest
    tests = functiontests(localfunctions);
end

function overallTest(testCase)
    global LEAF;
    del = 0.5;
    deldel = 1.5;
    for L = 5:15
        for n = 1:5
            state = RscaleSubTreeCoupled.unitTests.dummyState(ExpTree(L, 1));
            s = state.tree;
            for i = state.nodes
                nu = progeny(s, i, LEAF);
                [~, ~, t0Obs] = RscaleSubTreeCoupled.getNodePars(state, nu);
                [aObs, bObs] = RscaleSubTreeCoupled.getAgePars(...
                    state, i, t0Obs, del, deldel);

                t0Exp = min([s(nu(1, :)).time]);
                ti = s(i).time;
                aExp = t0Exp + del * (ti - t0Exp);
                bExp = t0Exp + (del + deldel) * (ti - t0Exp);

                assertEqual(testCase, t0Obs, t0Exp);
                assertEqual(testCase, aObs, aExp);
                assertEqual(testCase, bObs, bExp);
            end
        end
    end
end

function setupOnce(testCase)
    RscaleSubTreeCoupled.unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    RscaleSubTreeCoupled.unitTests.teardownOnce(testCase);
end

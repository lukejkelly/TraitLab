function tests = narrowOutputsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Compare new and old implementations
    for addClades = 0:1
        for h = 1:1e2
            L = 8 + ceil(rand * 8);
            state = dummyState(L, addClades);
            prior = struct('isclade', addClades);

            for iP = setdiff(state.nodes, state.root)
                for i = state.tree(iP).child
                    [jObs, iPObs, jPObs, OKObs] ...
                        = EchooseCoupled.narrowOutputs(i, state, prior);
                    assertEqual(testCase, iPObs, iP);
                    [jExp, iPExp, jPExp, OKExp] = narrowOutputsCheck(i, state, prior);

                    assertEqual(testCase, jObs, jExp);
                    assertEqual(testCase, iPObs, iPExp);
                    assertEqual(testCase, jPObs, jPExp);
                    assertEqual(testCase, OKObs, OKExp);
                end
            end
        end
    end
end

function [j, iP, jP, OK] = narrowOutputsCheck(i, state, prior)
    global OTHER
    s = state.tree;

    iP = s(i).parent;
    jP = s(iP).parent;
    j = s(jP).child(OTHER(s(iP).sibling));

    if (s(j).time < s(iP).time)
        if (~prior.isclade || (prior.isclade && isequal(s(iP).clade, s(jP).clade)))
            OK = 1;
        else
            OK = 0;
        end
    else
        OK = 0;
    end
end

% Dummy state
function state = dummyState(L, addClades)
    global ROOT
    state = unitTests.dummyState(ExpTree(L, 1e-3));
    if addClades
        s = state.tree;

        clade = synthclades(s, ceil(rand * L / 2), 2, 1 - rand^3);
        rootmax = 2 * s([s.type] == ROOT).time;
        prior = unitTests.clade2prior(clade, rootmax);

        state.tree = treeclades(state.tree, prior.clade);
        state = UpdateClades(state, [state.leaves, state.nodes], ...
                             size(prior.clade, 2));
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

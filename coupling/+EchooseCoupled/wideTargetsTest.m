function tests = wideTargetsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Compare new and old implementations
    for addClades = 0:1
        for i = 1:1e2
            L = 8 + ceil(rand * 8);
            state = dummyState(L, addClades);
            prior = struct('isclade', addClades);

            vObs = EchooseCoupled.wideTargets(state, prior);
            vExp = wideTargetsCheck(state, prior);

            assertFalse(testCase, any(isnan(vObs)));
            assertFalse(testCase, any(isnan(vExp)));
            assertEqual(testCase, sort(vObs), sort(vExp));
        end
    end
end

function v = wideTargetsCheck(state, prior)
    global ROOT OTHER
    Root = state.root;
    s = state.tree;
    N = 2 * state.NS - 1;

    u = nan(1, N);

    for i = 1:N
        if ( i==Root || ( s(s(i).parent).type==ROOT && s(s(s(i).parent).child(OTHER(s(i).sibling))).time<s(i).time ) ) % i=root or i=(root oldest child)
            u(i, :) = 0;
            continue;
        end
        for j = 1:N
            iP=s(i).parent;
            jP=s(j).parent;

            if i~=j && iP~=jP && i~=jP && j~=iP && s(j).time<s(iP).time && s(i).time<s(jP).time
                if (~prior.isclade ||(prior.isclade && isequal(s(iP).clade,s(jP).clade)))
                     u(i, j) = 1;
                else
                     u(i, j) = 0;
                end
            else
                u(i, j) = 0;
            end
        end
    end
    v = find(u);
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

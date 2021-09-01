function tests = GetLeavesTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    for L = 5:15
        state = dummyState(L);
        leaves = state.leaves;
        for m = 1:10
            leavesExp = leaves(randsample(L, ceil(rand * L)));
            names = {state.tree(leavesExp).Name};
            leavesObs = GetLeaves(state.tree, names);
            assertEqual(testCase, leavesObs, leavesExp);
        end
    end
end

function oldTest(testCase)
    for L = 5:15
        state = dummyState(L);
        for m = 1:10
            leaves = state.leaves(randsample(L, ceil(rand * L)));
            names = {state.tree(leaves).Name};
            leavesObs = GetLeaves(state.tree, names);
            leavesExp = oldGetLeaves(state.tree, names);
            assertEqual(testCase, leavesObs, leavesExp);
        end
    end
end

function y = oldGetLeaves(s, langs)
    GlobalSwitches;
    y=[];
    for k=1:size(langs,2)
        fi=find(strcmp(langs(k),{s.Name}));
        %if length(fi)~=1, disp('One of the Languages was not found or was not unique'); end
        y=[y,fi];
    end
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-4));
    state.rho = 3 / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

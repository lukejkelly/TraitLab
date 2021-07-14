function tests = rNeighbourTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global ROOT
    state = dummyState(8);
    s = state.tree;
    r = state.root;
    inds = find([s.type] < ROOT);
    n = 1e6;
    for i = inds
        [poss, q] = GetLegalCoupled.getPoss(s, i, r);
        new = arrayfun(@(~) GetLegalCoupled.rNeighbour(poss, q), (1:n)');
        p = mean(new == poss);
        hold on;
        plot(i, p - 1 / q, 'x');
        hold off;
    end
    xlabel('Branch');
    ylabel('Error in proportion of samples to each neighbour');
    fprintf('Error in observed proportions from %g samples\n', n);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for c = 1:(1 + poissrnd(3))
        state = AddCat(state);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 0;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

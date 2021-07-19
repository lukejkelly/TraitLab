function tests = sampleBranchByCatCountTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    n = 5e4;
    [new, old] = deal(nan(n, 1));
    state = dummyState(10);
    for i = 1:n
        new(i) = sampleBranchByCatCount(state);
        old(i) = sampleBranchByCatCountOld(state);
    end

    n_inds = 1:(2 * state.NS);
    p_new = mean(new == n_inds, 1);
    p_old = mean(old == n_inds, 1);

    p_exp = state.cat' ./ state.ncat;
    inds = find(p_exp);

    plot(n_inds(inds), [p_new(inds); p_old(inds)] - p_exp(inds), 'o');
    refline(0, mean((p_new(inds) - p_old(inds)).^2));
    xlabel('Node');
    ylabel('Difference from true proportion');
    legend('new', 'old', 'average squared difference');

    fprintf('Observed and expected proportions from %g samples\n', n);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);

end

function node = sampleBranchByCatCountOld(state)
    % Old function for drawing a branch proportional to its cat count
    r = ceil(state.ncat * rand);
    node = find(cumsum(state.cat) >= r, 1);
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

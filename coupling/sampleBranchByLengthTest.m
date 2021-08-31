function tests = sampleBranchByLengthTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    global ROOT;
    n = 5e4;
    [new, old] = deal(nan(n, 1));
    state = dummyState(10);
    for i = 1:n
        new(i) = sampleBranchByLength(state);
        old(i) = sampleBranchByLengthOld(state);
    end

    n_inds = 1:(2 * state.NS);
    p_new = mean(new == n_inds, 1);
    p_old = mean(old == n_inds, 1);

    p_exp = zeros(size(n_inds));
    s = state.tree;
    for i = n_inds
        if s(i).type < ROOT
            p_exp(i) = (s(s(i).parent).time - s(i).time) / state.length;
        end
    end

    plot(n_inds, [p_new; p_old] - p_exp, 'o');
    refline(0, mean((p_new - p_old).^2));
    xlabel('Node');
    ylabel('Difference from true proportion');
    legend('new', 'old', 'average squared difference');

    fprintf('Observed and expected proportions from %g samples\n', n);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function r = sampleBranchByLengthOld(state)
    % core/AddCat code for drawing a branch proportional to its length
    global ADAM ROOT
    L = 2 * state.NS;
    OK = 0;
    roottime = state.tree(state.root).time;
    while ~OK
        r = ceil(L * rand);
        if any(state.tree(r).type == [ADAM, ROOT])
            OK = 0;
        else
            dt = state.tree(state.tree(r).parent).time - state.tree(r).time;
            if rand < dt / roottime
                OK = 1;
            else
                OK = 0;
            end
        end
    end
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

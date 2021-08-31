function tests = GetLegalTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    rangeL = 6:10;
    n_i = length(rangeL);
    n_j = 1e2;
    for i = 1:n_i
        [old, new, q1, q2] = deal(zeros(1, n_j));
        state = dummyState(rangeL(i));
        s = state.tree;
        root = state.root;
        for j = 1:n_j
            old(j) = find(cumsum(state.cat) >= ceil(state.ncat * rand), 1);
            [new(j), q1(j), q2(j)] = GetLegal(s, old(j), root);
        end
        for u_old = unique(old)
            inds = find(old == u_old);
            assertTrue(testCase, all(q1(inds) == q1(inds(1))));
        end
        for u_new = unique(new)
            inds = find(new == u_new);
            assertTrue(testCase, all(q2(inds) == q2(inds(1))));
        end
    end
end

function distributionTest(testCase)
    global OTHER
    L = 10;
    n_i = 1e5;
    state = dummyState(L);
    s = state.tree;
    root = state.root;
    [cObs, cExp] = deal(zeros(1, 2 * L));
    for i = 1:n_i
        old = find(cumsum(state.cat) >= ceil(state.ncat * rand), 1);
        [new, ~, ~] = GetLegal(s, old, root);
        cObs(new) = cObs(new) + 1;
    end
    cObs = cObs / n_i;

    p = state.cat' ./ state.ncat;
    inds = find(p);
    for j = inds
        k = setdiff([s(j).parent, s(j).child, ...
                     s(s(j).parent).child(OTHER(s(j).sibling))], ...
                    root);
        cExp(k) = cExp(k) + p(j) / length(k);
    end

    fprintf('Proportion of matching samples after %g trials\n', n_i);
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs(inds)));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp(inds)));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs(inds) - cExp(inds)));
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for c = 1:(5 + poissrnd(3))
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
end

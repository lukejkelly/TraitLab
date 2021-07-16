function tests = sampleCatIndexTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Ensure that sampled index corresponds to a catastrophe
    M = 1e2;
    N = 1e2;
    for m = 1:M
        L = 5 + poissrnd(2);
        state = dummyState(L);
        k = nan(1, N);
        for n = 1:N
            [i, j] = sampleCatIndex(state);
            k(n) = state.tree(i).catloc(j);
        end
        assertLength(testCase, k, n);
    end
end

function distributionTest(testCase)
    % Ensure we are sampling uniformly across all the catastrophes on the tree
    M = 6;
    N = 1e5;
    for m = 1:M
        L = 5 + poissrnd(2);
        state = dummyState(L);
        [~, inds] = sampleCatIndexCoupled.getCats(state);
        fObs = zeros(1, state.ncat);
        for N = 1:N
            [i, j] = sampleCatIndex(state);
            k = find(i == inds(:, 1) & j == inds(:, 2));
            fObs(k) = fObs(k) + 1;
        end
        gObs = cumsum(fObs) / N;
        gExp = (1:state.ncat) / state.ncat;

        nexttile;
        plot(1:state.ncat, gObs, 'o', 1:state.ncat, gExp, 'x');
        legend('ECDF', 'CDF');
        xlabel('Catastrophe index');
        ylabel('(E)CDF');
        title(sprintf('x: %d leaves, %d cats', L, state.ncat));
    end
    fprintf('Observed and expected proportions from %g samples\n', N);
    v1 = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2 * rand));
    for c = 1:(3 + poissrnd(3))
        state = AddCat(state);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 1;
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

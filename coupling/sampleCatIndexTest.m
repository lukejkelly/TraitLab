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
            i = sampleBranchByCatCount(state);
            j = sampleCatIndex(state, i);
            k(n) = state.tree(i).catloc(j);
        end
        assertFalse(testCase, any(isnan(k)));
    end
end

function distributionTest(testCase)
    % Ensure for each branch we are sampling uniformly from its catastrophes
    M = 4;
    N = 1e4;
    for m = 1:M
        L = 5 + poissrnd(2);
        state = dummyState(L);
        [~, inds] = sampleCatIndex.getCats(state);
        C = unique(inds(:, 1)');
        for c = 1:length(C)
            nexttile;
            i = C(c);
            j = arrayfun(@(~) sampleCatIndex(state, i), 1:N);
            fObs = cumsum(mean((1:state.cat(i))' == j, 2));
            fExp = (1:state.cat(i)) / state.cat(i);
            plot(1:state.cat(i), fObs(:), 'o', 1:state.cat(i), fExp(:), 'x');
            xlabel(sprintf('i = %d : %d cats', i, state.cat(i)));
            ylabel('(E)CDF');
            title(sprintf('L = %d : %d cats', L, state.ncat));
        end
        legend('ECDF', 'CDF', 'location', 'southeast');
    end
    fprintf('Observed and expected proportions from %g samples\n', N);
    v1 = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2 * rand));
    for c = 1:(5 + poissrnd(3))
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

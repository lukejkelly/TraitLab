function tests = samplePriorCountsTest
    tests = functiontests(localfunctions);
end

function distributionTest(testCase)
    % Simulate catastrophes from prior, compare with relevant part of LogPrior
    % We check the support for combinations of 0 and 1 catastrophe
    global VARYRHO ANST

    [~, a, k] = LogRhoPrior(0);
    b = 1 / k;

    rangeL = 3:6;
    n_i = length(rangeL);
    n_j = 1e4;

    for i = 1:n_i
        L = rangeL(i);
        state = dummyState(L);

        p = 0:(2 * L);
        [c0, c1, q0, q1] = deal(zeros(length(p), 1));
        for j = 1:n_j
            VARYRHO = 0;
            cat0 = ResampleCatastrophes.samplePriorCounts(state);
            switch sum(cat0)
            case 0
                k = 1;
            case 1
                k = find(cat0) + 1;
            otherwise
                k = [];
            end
            c0(k) = c0(k) + 1;

            VARYRHO = 1;
            cat1 = ResampleCatastrophes.samplePriorCounts(state);
            switch sum(cat1)
            case 0
                k = 1;
            case 1
                k = find(cat1) + 1;
            otherwise
                k = [];
            end
            c1(k) = c1(k) + 1;
        end
        f0 = cumsum(c0) / n_j;
        f1 = cumsum(c1) / n_j;

        for j = [0, find([state.tree.type] <= ANST)]
            state.cat(:) = 0;
            if j > 0
                state.cat(j) = 1;
            end
            state.ncat = sum(state.cat);
            VARYRHO = 0;
            q0(j + 1) = exp(ResampleCatastrophes.logPrior(state));
            VARYRHO = 1;
            q1(j + 1) = exp(ResampleCatastrophes.logPrior(state));
        end
        h0 = cumsum(q0);
        h1 = cumsum(q1);

        subplot(n_i, 2, 2 * i - 1);
        plot(p, f0, 'o', p, h0, '*');
        title(sprintf('L  = %i: PP(%.3g)', L, state.rho));
        subplot(n_i, 2, 2 * i);
        plot(p, f1, 'o', p, h1, '*');
        title(sprintf('L  = %i: PG(%.3g, %.3g)', L, a, b / (b + state.length)));
        ylabel('CDF');
    end
    xlabel('Catastrophe sample index');
    legend('Sample', 'Expected');

    fprintf('Observed and expected proportions from %g samples\n', n_j);
    v1 = input('Are these distributions comparable? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);
end

% Dummy states
function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-4));
    state.rho = 2 / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT VARYRHO BORROWING
    BORROWING = 0;
    MCMCCAT = 1;
    testCase.TestData.VARYRHO = VARYRHO;
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    global VARYRHO
    VARYRHO = testCase.TestData.VARYRHO;
    close;
end

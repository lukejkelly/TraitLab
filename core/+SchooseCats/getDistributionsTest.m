function tests = getDistributionsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    L = 8;
    N = 1e4;
    state = dummyState(L);
    r = state.root;
    for j = 1:length(state.nodes)
        i = state.nodes(j);
        [~, kT, jT, a, b] = SchooseTime.nodeTimesAndRanges(i, state);
        if i == r
            newage = a + rand * (b - a);
        else
            newage = jT + rand * (kT - jT);
        end

        [catc, ncat, wc, wn] = SchooseCats.branchesAndWeights(state, i, newage);
        [rf, lf, lb] = SchooseCats.getDistributions(catc, ncat, wc, wn);

        x = nan(N, 3);
        y = nan(N, 1);
        for n = 1:N
            x(n, :) = rf();
            y(n) = ismembertol(lf(x(n, :)), log(mnpdf(x(n, :), wn)), 1e-10);
        end

        nexttile;
        plot(ncat * wn, mean(x, 1), 'o');
        refline(1, 0);
        title(sprintf('%i', i));
        xlabel('Expected proportion');
        ylabel('Observed proportion');

        assertTrue(testCase, all(y));
        assertEqual(testCase, lb, log(mnpdf(catc(:)', wc)));
    end
    v = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end

function setupOnce(testCase)
    global MCMCCAT
    unitTests.setupOnce(testCase);
    MCMCCAT = 1;
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for i = 1:(3 + poissrnd(3))
        state = AddCat(state);
    end
end

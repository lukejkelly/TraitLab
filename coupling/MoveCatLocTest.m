function tests = MoveCatLocTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    M = 1e3;
    [ncat, cat, i, U, OK, logq] = deal(nan(M, 1));
    for m = 1:M
        L = 8 + ceil(rand * 8);
        state = dummyState(L);
        [nstate, UObs, OKObs, logqObs] = MoveCatLoc(state);

        ncat(m) = nstate.ncat == state.ncat;
        cat(m) = all(nstate.cat == state.cat);
        if state.ncat == 0
            i(m) = true;
            U(m) = isempty(UObs);
            OK(m) = OKObs == 0;
        else
            [locs, inds] = sampleCatIndexCoupled.getCats(state);
            [locn, indn] = sampleCatIndexCoupled.getCats(nstate);

            ks = ~ismember(locs, locn);
            kn = ~ismember(locn, locs);
            i(m) = inds(ks, 1) == indn(kn, 1);

            U(m) = all(UObs == above(indn(kn, 1), nstate.tree, nstate.root));
            OK(m) = OKObs == 1;
        end
        logq(m) = logqObs == 0;
    end
    assertTrue(testCase, all(ncat));
    assertTrue(testCase, all(cat));
    assertTrue(testCase, all(i));
    assertTrue(testCase, all(U));
    assertTrue(testCase, all(OK));
    assertTrue(testCase, all(logq));
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, rand * 1e-2));
    for c = 1:poissrnd(2)
        state = AddCat(state);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

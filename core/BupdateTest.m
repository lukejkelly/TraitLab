function tests = BupdateTest
    tests = functiontests(localfunctions);
end

function catastropheTest(testCase)
    % mostly checking catastrophe-related outputs
    global BORROWING MCMCCAT NARROW WIDE
    M = 10;
    N = 10;
    for BORROWING = [0, 1]
        for MCMCCAT = [0, 1]
            for mt = [NARROW, WIDE]
                for m = 1:M
                    L = 5 + ceil(rand * 5);
                    [state, theta, prior] = dummyState(L);
                    for n = 1:N
                        newage = [];
                        while isempty(newage)
                            [i, j, k, newage, ~, ncat, cat, loc] = Bchoose(...
                                state, mt, theta, prior);
                        end
                        [p, q, h] = Bcats.pqh(state, i);
                        % j's parent is now k, h's [sib(i)] parent now pa(pa(i))
                        [nstate, ~, ~] = Bupdate(...
                            state, i, j, k, newage, ncat, cat, loc);

                        assertEqual(testCase, [nstate.tree([j, h, p]).parent], ...
                                    [p, q, k]);
                        if MCMCCAT
                            assertEqual(testCase, nstate.ncat, ncat);
                            assertEqual(testCase, nstate.cat([j, h, p]), ...
                                        structfun(@(x) x, cat));
                            if BORROWING
                                assertEqual(testCase, ...
                                            {nstate.tree([j, h, p]).catloc}, ...
                                            {sort(loc.j(1:cat.j)), ...
                                             sort(loc.h(1:cat.h)), ...
                                             sort(loc.p(1:cat.p))});
                            else
                                assertEqual(testCase, ...
                                            {nstate.tree([j, h, p]).catloc}, ...
                                            {[], [], []});
                            end
                        else
                            assertEqual(testCase, nstate.ncat, state.ncat);
                            assertEqual(testCase, nstate.cat, state.cat);
                            assertEqual(testCase, {nstate.tree.catloc}, ...
                                        {state.tree.catloc});
                        end

                    end
                end
            end
        end
    end
end

function [state, theta, prior] = dummyState(L)
    global ROOT
    theta = rand * 1e-2;

    s = ExpTree(L, theta);
    clade = synthclades(s, ceil(rand * L / 2), 2, 1 - rand^3);
    rootmax = (1 + rand) * s([s.type] == ROOT).time;
    prior = unitTests.clade2prior(clade, rootmax);

    state = unitTests.dummyState(s);
    state = UpdateClades(state, [state.leaves, state.nodes], ...
                         size(prior.clade, 2));
    for c = 1:poissrnd(4)
        state = AddCat(state);
    end
    state.rho = state.ncat / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    global VARYRHO
    unitTests.setupOnce(testCase);
    testCase.TestData.VARYRHO = VARYRHO;
    VARYRHO = 1;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

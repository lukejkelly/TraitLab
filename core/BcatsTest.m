function tests = BcatsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % check outputs, we have already checked Bcats.getDistributions for VARYRHO
    global BORROWING NARROW WIDE
    M = 20;
    N = 1e2;
    for BORROWING = [0, 1]
        for m = 1:M
            for addClades = 0:1
                L = 5 + ceil(rand * 5);
                [state, theta, prior] = dummyState(L, addClades);
                s = state.tree;
                r = state.root;
                for n = 1:N
                    if rand < 0.5
                        mt =  NARROW;
                    else
                        mt = WIDE;
                    end
                    newage = [];
                    while isempty(newage)
                        [i, j, k, newage, ~, ~, ~, ~] = Bchoose(...
                            state, mt, theta, prior);
                    end
                    [ncatObs, cat, loc, logqObs] = Bcats(state, i, j, k, newage);

                    [p, q, h] = Bcats.pqh(state, i);

                    ncatExp = state.ncat + sum(structfun(@(x) x, cat)) ...
                                  - sum(state.cat([j, h, p]));
                    assertEqual(testCase, ncatObs, ncatExp, 'AbsTol', 1e-12);

                    if BORROWING
                        assertEqual(testCase, structfun(@(x) x, cat), ...
                                    structfun(@length, loc));
                    else
                        assertNotEmpty(testCase, cat);
                        assertEmpty(testCase, loc);
                    end

                    dj = newage - s(j).time;
                    wj = dj / (s(k).time - s(j).time);
                    dh = s(p).time - s(h).time;
                    wh = dh / (s(q).time - s(h).time);
                    jc = state.cat(j);
                    hc = state.cat(h);
                    pc = state.cat(p);
                    switch r
                    case p
                        logqExp = Bcats.logPriorCounts(hc, dh, []) ...
                            - log(binopdf(cat.j, jc, wj));
                    case j
                        logqExp = log(binopdf(hc, hc + pc, wh)) ...
                            - Bcats.logPriorCounts(cat.j, dj, []);
                    otherwise
                        logqExp = log(binopdf(hc, hc + pc, wh)) ...
                            - log(binopdf(cat.j, jc, wj));
                    end
                    assertEqual(testCase, logqObs, logqExp, 'AbsTol', 1e-12);
                end
            end
        end
    end
end

function [state, theta, prior] = dummyState(L, addClades)
    global ROOT
    theta = rand * 1e-2;

    s = ExpTree(L, theta);
    if addClades
        numClades = ceil(rand * L / 2);
    else
        numClades = 0;
    end
    clade = synthclades(s, numClades, 2, 1 - rand^3);
    rootmax = (1 + rand) * s([s.type] == ROOT).time;
    prior = unitTests.clade2prior(clade, rootmax);

    state = unitTests.dummyState(s);
    state = UpdateClades(state, [state.leaves, state.nodes], ...
                         size(prior.clade, 2));
    for c = 1:poissrnd(4)
        state = AddCat(state);
    end
    state.rho = (0.5 + state.ncat) / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    global MCMCCAT VARYRHO
    unitTests.setupOnce(testCase);
    testCase.TestData.VARYRHO = VARYRHO;
    VARYRHO = 1;
    MCMCCAT = 1;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

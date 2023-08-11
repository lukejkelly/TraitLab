function tests = getDistributionsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % compare independent(-ish) implementations
    global NARROW WIDE VARYRHO
    M = 50;
    N = 20;
    O = 20;
    for m = 1:M
        L = 8 + ceil(rand * 8);
        for addClades = 0:1
            [state, theta, prior] = dummyState(L, addClades);
            s = state.tree;
            for mt = [NARROW, WIDE]
                for VARYRHO = [0, 1]
                    for n = 1:N
                        newage = [];
                        while isempty(newage)
                            [i, j, k, newage, ~] = Bchoose(state, mt, theta, prior);
                        end
                        [rfObs, lfObs, lbObs] = Bcats.getDistributions(...
                            state, i, j, k, newage);

                        [~, a, bInv] = LogRhoPrior(0);
                        b = 1 / bInv;
                        [p, q, h] = Bcats.pqh(state, i);

                        jc = state.cat(j);
                        dj = newage - s(j).time;
                        if j == state.root
                            if VARYRHO
                                rfExp = @() PoissonGammaSample(a, b / (b + dj));
                                lfExp = @(x) PoissonGammaLogProb(x, a, b / (b + dj));
                            else
                                rfExp = @() poissrnd(state.rho * dj);
                                lfExp = @(x) log(poisspdf(x, state.rho * dj));
                            end
                        else
                            wj = dj / (s(k).time - s(j).time);
                            rfExp = @() binornd(jc, wj);
                            lfExp = @(x) log(binopdf(x, jc, wj));
                        end

                        [oObs, oExp] = deal(nan(O, 2));
                        for o = 1:O
                            tObs = rfObs();
                            tExp = rfExp();
                            oObs(o, :) = [lfObs(tObs), lfObs(tExp)];
                            oExp(o, :) = [lfExp(tObs), lfExp(tExp)];
                        end
                        assertEqual(testCase, oObs, oExp, 'AbsTol', 1e-12);

                        dh = s(p).time - s(h).time;
                        hc = state.cat(h);
                        if p == state.root
                            if VARYRHO
                                lbExp = PoissonGammaLogProb(hc, a, b / (b + dh));
                            else
                                lbExp = log(poisspdf(hc, dh * state.rho));
                            end
                        else
                            wh = dh / (s(q).time - s(h).time);
                            lbExp = log(binopdf(hc, hc + state.cat(p), wh));
                        end
                        assertEqual(testCase, lbObs, lbExp, 'AbsTol', 1e-12);
                    end
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
    for c = 1:poissrnd(3)
        state = AddCat(state);
    end
    state.rho = (0.5 + state.ncat) / state.length;
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING VARYRHO
    testCase.TestData.VARYRHO = VARYRHO;
    MCMCCAT = 1;
    BORROWING = 1;
end

function teardownOnce(testCase)
    global VARYRHO
    unitTests.teardownOnce(testCase);
    VARYRHO = testCase.TestData.VARYRHO;
end

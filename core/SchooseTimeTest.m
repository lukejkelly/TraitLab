function tests = SchooseTimeTest
    tests = functiontests(localfunctions);
end

function parameterTest(testCase)
    for L = 6:10
        state = dummyState(L);
        visited = zeros(1, state.NS - 1);

        while all(visited < 1e2)
            [i, newage, logq] = SchooseTime(state);

            visited(state.nodes == i) = visited(state.nodes == i) + 1;
            iT = state.tree(i).time;
            jT = max([state.tree(state.tree(i).child).time]);
            if i == state.root
                assertTrue(testCase, ...
                           mean([iT, jT]) <= newage && newage <= 2 * iT - jT);
                assertTrue(testCase, -log(2) <= logq && logq <= log(2));
                assertEqual(testCase, ...
                            logq, log(iT - jT) - log(newage - jT), ...
                            'AbsTol', 1e-12)
            else
                kT = state.tree(state.tree(i).parent).time;
                assertTrue(testCase, jT <= iT && iT <= kT);
                assertTrue(testCase, jT <= newage && newage <= kT);
                assertEqual(testCase, logq, 0);
            end
        end
    end
end

function marginalTest(testCase)
    L = 8;
    state = dummyState(L);
    r = state.root;

    N = 1e3;
    [iObs, iExp, newageObs, newageExp, logqObs, logqExp] = deal(zeros(N, 1));

    for n = 1:N
        [iObs(n), newageObs(n), logqObs(n)] = SchooseTime(state);
        [iExp(n), newageExp(n), logqExp(n)] = oldSchooseTime(state);
    end

    assertTrue(testCase, all(logqObs(iObs ~= r) == 0));
    assertTrue(testCase, all(logqExp(iExp ~= r) == 0));

    clf;
    for j = 1:length(state.nodes)
        i = state.nodes(j);
        nexttile;
        [nObs, eObs] = ecdf(newageObs(iObs == i));
        [nExp, eExp] = ecdf(newageExp(iExp == i));
        plot(eObs, nObs, eExp, nExp);
        xlabel('newage');
        title(sprintf('%i', i))
    end
    nexttile;
    [nObs, eObs] = ecdf(logqObs(iObs == r));
    [nExp, eExp] = ecdf(logqExp(iExp == r));
    plot(eObs, nObs, eExp, nExp);
    xlabel('logq');
    title(sprintf('%i', i))
    legend('Obs', 'Exp', 'Location', 'southeast');
    v = input('Do the CDFs in the figure match? Reply 1 for yes... ');
    assertEqual(testCase, v, 1);
end

function [i, newage, logq] = oldSchooseTime(state)
    % Old version of Schoose after removing comments and before reworking
    global ROOT

    s = state.tree;
    i = state.nodes(ceil(rand * (state.NS - 1)));

    iT = s(i).time;
    k = s(i).parent;
    j1 = s(i).child(1);
    j2 = s(i).child(2);
    jT = max(s(j1).time, s(j2).time);

    if s(i).type == ROOT
        tau = iT - jT;
        delta = (-0.5 + rand * 1.5) * tau;
        taup = tau + delta;
        logq = log(tau / taup);
        newage = jT + taup;
    else
        newage = jT + rand * (s(k).time - jT);
        logq = 0;
    end
end

function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
end

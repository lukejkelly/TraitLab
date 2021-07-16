function tests = EchooseCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global NARROW WIDE
    nR = 1e2;
    rangeL = 6:10;
    nL = length(rangeL);
    for mt = [NARROW, WIDE]
        fprintf('Move = %d\n', mt);
        for addClades = 0:1
            fprintf('\tClades = %d\n', addClades);
            prior = struct('isclade', addClades);
            [i_x, i_y, j_x, j_y, iP_x, iP_y, jP_x, jP_y, logq_x, logq_y, ...
                OK_x, OK_y] = deal(nan(nL, nR));
            for m = 1:nL
                L = rangeL(m);
                for n = 1:nR
                    [state_x, state_y] = coupledStates(L, addClades);
                    [i_x(m, n), i_y(m, n), j_x(m, n), j_y(m, n), ...
                     iP_x(m, n), iP_y(m, n), jP_x(m, n), jP_y(m, n), ...
                     logq_x(m, n), logq_y(m, n), OK_x(m, n), OK_y(m, n)] ...
                        = EchooseCoupled(state_x, state_y, mt, prior);
                end
            end
            assertFalse(testCase, any(isnan(i_x), 'all'));
            assertFalse(testCase, any(isnan(i_y), 'all'));
            assertFalse(testCase, any(isnan(j_x), 'all'));
            assertFalse(testCase, any(isnan(j_y), 'all'));
            assertFalse(testCase, any(isnan(iP_x), 'all'));
            assertFalse(testCase, any(isnan(iP_y), 'all'));
            assertFalse(testCase, any(isnan(jP_x), 'all'));
            assertFalse(testCase, any(isnan(jP_y), 'all'));
            assertEqual(testCase, logq_x, zeros(nL, nR));
            assertEqual(testCase, logq_y, zeros(nL, nR));
            assertFalse(testCase, any(isnan(OK_x), 'all'));
            assertFalse(testCase, any(isnan(OK_y), 'all'));

            assertTrue(testCase, any(OK_x, 'all'));
            assertTrue(testCase, any(OK_y, 'all'));

            assertEqual(testCase, i_x, i_y);
            assertEqual(testCase, j_x, j_y);
            assertEqual(testCase, iP_x, iP_y);
            assertEqual(testCase, jP_x, jP_y);
            assertEqual(testCase, OK_x, OK_y);
        end
    end
end

function marginalTest(testCase)
    global NARROW WIDE
    nR = 2e4;
    L = 10;
    for mt = [NARROW, WIDE]
        fprintf('Move = %d\n', mt);
        for addClades = 0:1
            fprintf('\tClades = %d\n', addClades);
            prior = struct('isclade', addClades);
            [state_x, state_y] = housekeptStates(L, addClades);

            [i_xC, j_xC, iP_xC, jP_xC, logq_xC, OK_xC] = deal(nan(nR, 1));
            [i_yC, j_yC, iP_yC, jP_yC, logq_yC, OK_yC] = deal(nan(nR, 1));

            [i_xM, j_xM, iP_xM, jP_xM, logq_xM, OK_xM] = deal(nan(nR, 1));
            [i_yM, j_yM, iP_yM, jP_yM, logq_yM, OK_yM] = deal(nan(nR, 1));

            for n = 1:nR
                [i_xC(n), i_yC(n), j_xC(n), j_yC(n), ...
                 iP_xC(n), iP_yC(n), jP_xC(n), jP_yC(n), ...
                 logq_xC(n), logq_yC(n), OK_xC(n), OK_yC(n)] ...
                    = EchooseCoupled(state_x, state_y, mt, prior);

                [i_xM(n), j_xM(n), iP_xM(n), jP_xM(n), logq_xM(n), OK_xM(n)] ...
                    = Echoose(state_x, mt, prior);
                [i_yM(n), j_yM(n), iP_yM(n), jP_yM(n), logq_yM(n), OK_yM(n)] ...
                    = Echoose(state_y, mt, prior);
            end

            assertEqual(testCase, logq_xC, zeros(nR, 1));
            assertEqual(testCase, logq_yC, zeros(nR, 1));
            assertEqual(testCase, logq_xM, zeros(nR, 1));
            assertEqual(testCase, logq_yM, zeros(nR, 1));

            assertTrue(testCase, any(OK_xC, 'all'));
            assertTrue(testCase, any(OK_yC, 'all'));
            assertTrue(testCase, any(OK_xM, 'all'));
            assertTrue(testCase, any(OK_yM, 'all'));

            % Filter by successful moves
            i_xC = i_xC(OK_xC == 1);
            j_xC = j_xC(OK_xC == 1);
            iP_xC = iP_xC(OK_xC == 1);
            jP_xC = jP_xC(OK_xC == 1);

            i_yC = i_yC(OK_yC == 1);
            j_yC = j_yC(OK_yC == 1);
            iP_yC = iP_yC(OK_yC == 1);
            jP_yC = jP_yC(OK_yC == 1);

            i_xM = i_xM(OK_xM == 1);
            j_xM = j_xM(OK_xM == 1);
            iP_xM = iP_xM(OK_xM == 1);
            jP_xM = jP_xM(OK_xM == 1);

            i_yM = i_yM(OK_yM == 1);
            j_yM = j_yM(OK_yM == 1);
            iP_yM = iP_yM(OK_yM == 1);
            jP_yM = jP_yM(OK_yM == 1);

            ui_x = unique(union(i_xC', i_xM'));
            for i_x = ui_x
                u_xC = unique(iP_xC(i_xC == i_x));
                u_xM = unique(iP_xM(i_xM == i_x));
                assertLength(testCase, u_xC, 1);
                assertEqual(testCase, u_xC, u_xM);
            end

            ui_y = unique(union(i_yC', i_yM'));
            for i_y = ui_y
                u_yC = unique(iP_yC(i_yC == i_y));
                u_yM = unique(iP_yM(i_yM == i_y));
                assertLength(testCase, u_yC, 1);
                assertEqual(testCase, u_yC, u_yM);
            end

            uj_x = unique(union(j_xC', j_xM'));
            for j_x = uj_x
                u_xC = unique(jP_xC(j_xC == j_x));
                u_xM = unique(jP_xM(j_xM == j_x));
                assertLength(testCase, u_xC, 1);
                assertEqual(testCase, u_xC, u_xM);
            end

            uj_y = unique(union(j_yC', j_yM'));
            for j_y = uj_y
                u_yC = unique(jP_yC(j_yC == j_y));
                u_yM = unique(jP_yM(j_yM == j_y));
                assertLength(testCase, u_yC, 1);
                assertEqual(testCase, u_yC, u_yM);
            end

            nexttile;
            plot(ui_x, cumsum(mean(i_xC == ui_x)), '*', ...
                 ui_x, cumsum(mean(i_xM == ui_x)), 'o')
            title(sprintf('m %d / c %d / $ i_x $', mt, addClades));

            nexttile;
            plot(uj_x, cumsum(mean(j_xC == uj_x)), '*', ...
                 uj_x, cumsum(mean(j_xM == uj_x)), 'o')
            title(sprintf('m %d / c %d / $ j_x $', mt, addClades));

            nexttile;
            plot(ui_y, cumsum(mean(i_yC == ui_y)), '*', ...
                 ui_y, cumsum(mean(i_yM == ui_y)), 'o')
            title(sprintf('m %d / c %d / $ i_y $', mt, addClades));

            nexttile;
            plot(uj_y, cumsum(mean(j_yC == uj_y)), '*', ...
                 uj_y, cumsum(mean(j_yM == uj_y)), 'o')
            title(sprintf('m %d / c %d / $ j_y $', mt, addClades));
        end
    end
    fprintf('Observed and expected proportions from %g samples\n', nR);
    v = input('Are these plots okay? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

% Dummy states
function [state_x, state_y] = coupledStates(L, addClades)
    global ROOT
    [state_x, state_y] = unitTests.coupledStates(L, 1e-3);
    if addClades
        s = state_x.tree;

        clade = synthclades(s, ceil(rand * L / 2), 2, 1 - rand^3);
        rootmax = 2 * s([s.type] == ROOT).time;
        prior = unitTests.clade2prior(clade, rootmax);

        state_x.tree = treeclades(state_x.tree, prior.clade);
        state_x = UpdateClades(state_x, [state_x.leaves, state_x.nodes], ...
                               size(prior.clade, 2));

        state_y.tree = treeclades(state_y.tree, prior.clade);
        state_y = UpdateClades(state_y, [state_y.leaves, state_y.nodes], ...
                               size(prior.clade, 2));
    end
end

function [state_x, state_y] = housekeptStates(L, addClades)
    global ROOT
    if addClades

        state = unitTests.dummyState(ExpTree(L, 1e-3));
        s = state.tree;
        clade = synthclades(s, ceil(rand * L / 2), 2, 1 - rand^3);
        rootmax = (1 + rand) * s([s.type] == ROOT).time;
        prior = unitTests.clade2prior(clade, rootmax);

        s_x = RandCladeTree(clade, {s(state.leaves).Name}, rootmax, 0);
        s_x = treeclades(s_x, prior.clade);
        state_x = unitTests.dummyState(s_x);
        state_x = UpdateClades(state_x, [state_x.leaves, state_x.nodes], ...
                               size(prior.clade, 2));

        s_y = RandCladeTree(clade, {s(state.leaves).Name}, rootmax, 0);
        s_y = treeclades(s_y, prior.clade);
        state_y = unitTests.dummyState(s_y);
        state_y = UpdateClades(state_y, [state_y.leaves, state_y.nodes], ...
                               size(prior.clade, 2));

        state_y = housekeeping(state_x, state_y);

    else
       [state_x, state_y] = unitTests.housekeptStates(L, 5e-3 * rand);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    global BORROWING MCMCCAT
    unitTests.setupOnce(testCase);
    BORROWING = 0;
    MCMCCAT = 0;
    clf;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

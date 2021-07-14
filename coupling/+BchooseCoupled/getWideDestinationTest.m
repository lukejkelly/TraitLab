function tests = getWideDestinationTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global ROOT
    for L = 6:10
        for addClades = 0:1
            [state_x, state_y] = coupledStates(L, addClades);
            s_x = state_x.tree;
            s_y = state_y.tree;
            for i = find([s_x.type] < ROOT)
                if addClades
                    r = BchooseCoupled.getWideCandidatesClade(i, s_x);
                    N = length(r);
                else
                    r = find([s_x.type] <= ROOT);
                    N = length(r);
                end

                if N > 4
                    [j_x, j_y, k_x, k_y, FAIL_x, FAIL_y] ...
                        = BchooseCoupled.getWideDestination(i, r, s_x, s_y);
                end
                assertEqual(testCase, j_x, j_y);
                assertEqual(testCase, k_x, k_y);
                assertEqual(testCase, FAIL_x, FAIL_y);
            end
        end
    end
end

function couplingTest(testCase)
    global ROOT
    nSamp = 1e4;
    rangeL = 8:2:12;
    for L = rangeL
        for addClades = 0:1
            [cObs, pExp] = deal(zeros(1, 2 * L - 1));
            [state_x, state_y] = housekeptStates(L, addClades);
            s_x = state_x.tree;
            s_y = state_y.tree;
            for i = find([s_x.type] < ROOT)
                if addClades
                    r = BchooseCoupled.getWideCandidatesClade(i, s_x);
                    N = length(r);
                else
                    r = find([s_x.type] <= ROOT);
                    N = length(r);
                end

                if N > 4
                    for nS = 1:nSamp
                        [j_x, j_y, ~, ~, ~, ~] ...
                            = BchooseCoupled.getWideDestination(i, r, s_x, s_y);
                        if j_x == j_y
                            cObs(i) = cObs(i) + 1;
                        end
                    end
                    z_x = BchooseCoupled.getWideDestination.valid(i, r, s_x);
                    z_y = BchooseCoupled.getWideDestination.valid(i, r, s_y);
                    pExp(i) = length(intersect(z_x, z_y)) / max(length(z_x), length(z_y));
                end
            end
            pObs = cObs / nSamp;
            nexttile;
            plot(1:(2 * L - 1), pObs, 'o', 1:(2 * L - 1), pExp, 'x');
            legend('Observed', 'Expected');
            xlabel('Branch i');
            ylabel('Proportion j_x = j_y');
            title(sprintf('L = %d : clades = %d', L, addClades));
        end
    end
    fprintf('Observed and expected proportions from %g samples\n', nSamp);
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

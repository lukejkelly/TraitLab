function tests = ResampleCatastrophesCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global BORROWING VARYRHO
    for L = 5:10
        for i = 1:10
            for BORROWING = 0:1
                for VARYRHO = 0:1
                    [state_x, state_y] = coupledStates(L);
                    [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
                            = ResampleCatastrophesCoupled(state_x, state_y);
                    assertTrue(testCase, ...
                               equaltrees(nstate_x.tree, nstate_y.tree) == 1);
                    assertEqual(testCase, nstate_x.cat, nstate_y.cat);
                    assertEqual(testCase, nstate_x.ncat, nstate_y.ncat);
                    if BORROWING
                        for j = 1:length(nstate_x.tree)
                            assertEqual(testCase, nstate_x.tree(j).catloc, ...
                                        nstate_y.tree(j).catloc);
                        end
                    end
                    assertEqual(testCase, U_x, U_y);
                    assertEqual(testCase, OK_x, OK_y);
                    assertEqual(testCase, logq_x, logq_y);
                end
            end
        end
    end
end

function overlapTest(testCase)
    % Compare  distributions restricted to 0 or 1 catastrophe
    global VARYRHO BORROWING
    BORROWING = 0;
    L = 10;
    n_i = 6;
    n_j = 1e4;
    [cObs, cExp] = deal(nan(n_i, 2));
    for VARYRHO = 0:1
        for i = 1:n_i
            [state_x, state_y] = housekeptStates(L);
            c = 0;
            for j = 1:n_j
                [nstate_x, nstate_y] = ResampleCatastrophesCoupled(state_x, ...
                                                                   state_y);
                if nstate_x.ncat <= 1 && all(nstate_x.cat == nstate_y.cat)
                    c = c + 1;
                end
            end
            cObs(i, VARYRHO + 1) = c / n_j;

            % evaluate overlap up to largest mean + 3 standard deviations
            dx = getBranchLengths(nstate_x);
            dy = getBranchLengths(nstate_y);
            Dx = nstate_x.length;
            Dy = nstate_y.length;
            if VARYRHO
                [~, a, k] = LogRhoPrior(0);
                b = 1 / k;

                px = b / (b + Dx);
                py = b / (b + Dy);

                qx = dx / (b + Dx);
                qy = dy / (b + Dy);

                ol = sum(min(px^a * [1, a * qx], py^a * [1, a * qy]));
            else
                zx = state_x.rho * dx;
                zy = state_y.rho * dy;
                Zx = state_x.rho * Dx;
                Zy = state_y.rho * Dy;

                ol = sum(min(exp(-Zx) * [1, zx], exp(-Zy) * [1, zy]));
            end
            cExp(i, VARYRHO + 1) = ol;
        end
    end

    fprintf('Proportion of matching samples in each of %g trials\n', n_j);
    fprintf('Repeated for %g 10-leaf trees\n\n', n_i);

    fprintf('VARYRHO = 0 : Poisson counts\n')
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs(:, 1)));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp(:, 1)));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs(:, 1) - cExp(:, 1)));
    v1 = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v1 == 1);

    fprintf('VARYRHO = 1 : Poisson-Gamma counts\n')
    fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs(:, 2)));
    fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp(:, 2)));
    fprintf('Difference: %s\n', sprintf('%+-3.2e\t', cObs(:, 2) - cExp(:, 2)));
    v2 = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v2 == 1);
end

function locationTest(testCase)
    % Check catastrophe locations coupled as expected
    global BORROWING VARYRHO
    BORROWING = 1;
    n_i = 10;
    n_j = 10;
    for VARYRHO = 0:1
        for i = 1:n_i
            [state_x, state_y] = housekeptStates(2 + i);
            for j = 1:n_j
                [nstate_x, nstate_y] ...
                    = ResampleCatastrophesCoupled(state_x, state_y);
                assertEqual(testCase, nstate_x.cat(:)', ...
                            cellfun(@length, {nstate_x.tree.catloc}));
                assertEqual(testCase, nstate_y.cat(:)', ...
                            cellfun(@length, {nstate_y.tree.catloc}));
                assertEqual(testCase, ...
                            min(nstate_x.cat(:)', nstate_y.cat(:)'), ...
                            cellfun(@(x, y) length(intersect(x, y)), ...
                                    {nstate_x.tree.catloc}, ...
                                    {nstate_y.tree.catloc}));
                % for k = 1:length(nstate_x.tree)
                %     cx = nstate_x.cat(k);
                %     lx = nstate_x.tree(k).catloc;
                %     assertEqual(testCase, cx, length(lx));
                %     cy = nstate_y.cat(k);
                %     ly = nstate_y.tree(k).catloc;
                %     assertEqual(testCase, cy, length(ly));
                %     assertEqual(testCase, min(cx, cy), length(intersect(lx, ly)));
                % end
            end
        end
    end
end

% Dummy states
function [state_x, state_y] = housekeptStates(L)
    [state_x, state_y] = unitTests.housekeptStates(L, 5e-3 * rand);
    state_x.rho = 3 * rand / state_x.length;
    state_y.rho = 3 * rand / state_y.length;
    for c = 1:(3 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
    for c = 1:(3 + poissrnd(3))
        state_x = AddCat(state_x);
    end
    for c = 1:(3 + poissrnd(3))
        state_y = AddCat(state_y);
    end
end

function [state_x, state_y] = coupledStates(L)
    [state_x, state_y] = unitTests.coupledStates(L, 5e-3 * rand);
    state_x.rho = 3 / state_x.length;
    state_y.rho = 3 / state_y.length;
    for c = 1:(3 + poissrnd(3))
        [state_x, state_y] = AddCatCoupled(state_x, state_y);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global BORROWING MCMCCAT VARYRHO
    BORROWING = 0;
    MCMCCAT = 1;
    testCase.TestData.VARYRHO = VARYRHO;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    global VARYRHO
    VARYRHO = testCase.TestData.VARYRHO;
end

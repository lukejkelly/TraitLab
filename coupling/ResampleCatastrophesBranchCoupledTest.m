function tests = ResampleCatastrophesBranchCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global BORROWING VARYRHO
    for BORROWING = 0:1
        for VARYRHO = 0:1
            for L = 5:10
                for i = 1:10
                    [state_x, state_y] = coupledStates(L);
                    [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
                            = ResampleCatastrophesBranchCoupled(state_x, state_y);
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
    M = 6;
    N = 1e4;
    [cObs, cExp] = deal(nan(M, 2));
    for VARYRHO = 0:1
        for m = 1:M
            [state_x, state_y] = housekeptStates(L);
            c = 0;
            for n = 1:N
                [nstate_x, nstate_y, U_x, U_y] = ResampleCatastrophesBranchCoupled(...
                    state_x, state_y);
                if U_x(1) == U_y(1) && nstate_x.cat(U_x(1)) == nstate_y.cat(U_y(1))
                    c = c + 1;
                end
            end
            cObs(m, VARYRHO + 1) = c / N;

            % evaluate overlap up to largest mean + 3 standard deviations
            dx = getBranchLengths(nstate_x);
            dy = getBranchLengths(nstate_y);
            Dx = nstate_x.length;
            Dy = nstate_y.length;

            % probabilty of picking same branch under maximal coupling
            p_i = min(dx / Dx, dy / Dy);

            % probability of sampling same count given same branch
            f_i = zeros(size(p_i));
            for i = find(p_i)
                f_i(i) = sum(exp(arrayfun(...
                    @(z) min(ResampleCatastrophesBranch.logPriorCount(state_x, i, z), ...
                             ResampleCatastrophesBranch.logPriorCount(state_y, i, z)), ...
                    0:(2 * max(nstate_x.ncat, nstate_y.ncat)))));
            end
            cExp(m, VARYRHO + 1) = sum(p_i .* f_i);
        end
    end

    fprintf('Proportion of matching samples in %g trials\n', N);
    fprintf('Repeated for %g 10-leaf trees\n\n', M);

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
                [nstate_x, nstate_y, U_x, U_y] ...
                    = ResampleCatastrophesBranchCoupled(state_x, state_y);
                i_x = U_x(1);
                i_y = U_y(1);
                assertEqual(testCase, nstate_x.cat(i_x), ...
                            length(nstate_x.tree(i_x).catloc));
                assertEqual(testCase, nstate_y.cat(i_y), ...
                            length(nstate_y.tree(i_y).catloc));
                assertEqual(testCase, ...
                            min(nstate_x.cat(i_x), nstate_y.cat(i_y)), ...
                            length(intersect(nstate_x.tree(i_x).catloc, ...
                                             nstate_y.tree(i_y).catloc)));
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
    global MCMCCAT VARYRHO
    MCMCCAT = 1;
    testCase.TestData.VARYRHO = VARYRHO;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    global VARYRHO
    VARYRHO = testCase.TestData.VARYRHO;
end

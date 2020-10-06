function tests = patternCountsTest
    % Unit-testing likelihood calculations in patternCounts
    tests = functiontests(localfunctions);
end

%% Test Functions
function Data3Missing0Registration0Test(testCase)
    testState = load('borrowing/patternCountsTestData3Missing0Registration0');
    Missing0Registration0(testCase, testState);
end

function Data3Missing0Registration1Test(testCase)
    testState = load('borrowing/patternCountsTestData3Missing0Registration1');
    Missing0Registration1(testCase, testState);
end

function Data3Missing1Registration0Test(testCase)
    testState = load('borrowing/patternCountsTestData3Missing1Registration0');
    Missing1Registration0(testCase, testState);
end

% Helper functions
function Missing0Registration0(testCase, testVars)
    % Although we are not modelling missing data, traits with missing entries do
    % not get removed, other than those which are potentially absent in every
    % taxon, but are accounted for by observedPatternCounts and patternCounts

    % Test domain
    global LOSTONES MISDAT;

    % Loading test state
    state = testVars.state;
    MISDAT = testVars.MISDAT;
    LOSTONES = testVars.LOSTONES;

    % Comparing test domain and loaded state
    assertEqual(testCase, MISDAT, 0);
    assertEqual(testCase, LOSTONES, 0);

    % Initialising other parameters
    borPars = borrowingParameters(state.NS);
    rl = state.leaves;
    eta =  2^state.NS - 1;
    x_T = gamrnd(1, 1, eta, 1);

    % No registration or missing data correction to sum of pattern frequencies
    C = 0;

    % Checking output for uniform data
    Missing0LogLkd2Check(testCase, state, rl, x_T, borPars, ones(eta, 1), C);

    % Checking output for non-uniform data
    indReps = poissrnd(2, 1, state.L);
    state.L = sum(indReps);
    for i = state.leaves
        state.tree(i).dat = repelem(state.tree(i).dat, indReps);
    end
    dat = reshape([state.tree(fliplr(rl)).dat], state.L, state.NS);
    dat(any(dat == 2, 2), :) = [];
    n_T = sum(bi2de(dat, 'left-msb') == (1:eta), 1)';
    Missing0LogLkd2Check(testCase, state, rl, x_T, borPars, n_T, C);
end

function Missing0Registration1(testCase, testVars)
    % Traits attested in 1 leaf or fewer removed

    % Test domain
    global LOSTONES MISDAT;

    % Loading test state
    state = testVars.state;
    MISDAT = testVars.MISDAT;
    LOSTONES = testVars.LOSTONES;

    % Comparing test domain and loaded state
    assertEqual(testCase, MISDAT, 0);
    assertEqual(testCase, LOSTONES, 1);

    % Initialising other parameters
    borPars = borrowingParameters(state.NS);
    rl = state.leaves;
    eta =  2^state.NS - 1;
    x_T = gamrnd(1, 1, eta, 1);

    % Registration correction to sum of pattern frequencies
    indsReg = borPars(end).S > 1;
    C = sum(x_T(~indsReg));

    % Checking output for uniform data
    Missing0LogLkd2Check(testCase, state, rl, x_T, borPars, double(indsReg), C);

    % Checking output for non-uniform data
    indReps = poissrnd(2, 1, state.L);
    state.L = sum(indReps);
    for i = state.leaves
        state.tree(i).dat = repelem(state.tree(i).dat, indReps);
    end
    dat = reshape([state.tree(fliplr(rl)).dat], state.L, state.NS);
    dat(any(dat == 2, 2), :) = [];
    n_T = sum(bi2de(dat, 'left-msb') == (1:eta), 1)';
    Missing0LogLkd2Check(testCase, state, rl, x_T, borPars, n_T, C);
end

function Missing0LogLkd2Check(testCase, state, rl, x_T, borPars, n_T, C)
    % x_T and n_T were computed on a tree with ordering rl and indexed by binary
    % representation in P_b
    clear patternCounts;
    [intLogLikObs, logLikObs] = patternCounts(state, rl, x_T, borPars);
    intLogLikExp = sum(n_T .* log(x_T)) - sum(n_T) * log(sum(x_T) - C);
    logLikExp = sum(n_T .* log(x_T * state.lambda)) ...
                - state.lambda * (sum(x_T) - C);
    assertEqual(testCase, intLogLikObs, intLogLikExp, 'AbsTol', 1e-10);
    assertEqual(testCase, logLikObs, logLikExp, 'AbsTol', 1e-10);
end

function Missing1Registration0(testCase, testVars)
    % Only traits which are not present in at least one leaf get discarded

    % Test domain
    global LOSTONES MISDAT;

    % Loading test state
    state = testVars.state;
    MISDAT = testVars.MISDAT;
    LOSTONES = testVars.LOSTONES;

    % Comparing test domain and loaded state
    assertEqual(testCase, MISDAT, 1);
    assertEqual(testCase, LOSTONES, 0);

    % Initialising other parameters
    borPars = borrowingParameters(state.NS);
    rl = state.leaves;
    eta =  2^state.NS - 1;
    x_T = gamrnd(1, 1, eta, 1);

    % Correction for patterns not observed at at least one leaf
    xi = fliplr([state.tree(rl).xi]);
    C = sum(x_T .* prod(bsxfun(@power, 1 - xi, borPars(end).P_b), 2));

    % Checking output for uniform data
    n_U = zeros(3^state.NS - 1, 1);
    n_U(reshape([state.tree(fliplr(rl)).dat], [], 3) ...
       * 3.^((state.NS - 1):(-1):0)') = 1;
    Missing1LogLkd2Check(testCase, state, rl, x_T, borPars, n_U, C);

    % Checking output for non-uniform data
    indReps = poissrnd(2, 1, state.L);
    state.L = sum(indReps);
    for i = state.leaves
        state.tree(i).dat = repelem(state.tree(i).dat, indReps);
    end
    dat = reshape([state.tree(fliplr(rl)).dat], state.L, state.NS);
    n_T = sum(dat * 3.^((state.NS - 1):(-1):0)' == (1:(3^state.NS - 1)), 1)';
    Missing1LogLkd2Check(testCase, state, rl, x_T, borPars, n_T, C);
end

function Missing1LogLkd2Check(testCase, state, rl, x_T, borPars, n_T, C)
    % x_T and n_T were computed on a tree with ordering rl
    % x_T indexed by binary representation in P_b
    % n_T indexed by ternary representation
    clear patternCounts;
    [intLogLikObs, logLikObs] = patternCounts(state, rl, x_T, borPars);
    xi = fliplr([state.tree(rl).xi]);
    % Expanded pattern set
    Q = zeros(state.NS^3 - 1, state.NS);
    for j = 1:state.NS
        Qj = repmat(repelem(0:2, 3^(state.NS - j))', 3^(j - 1), 1);
        Q(:, j) = Qj(2:end);
    end
    [intLogLikExp_i, logLikExp_i, y_T] = deal(zeros(size(n_T)));
    for i = 1:length(n_T)
        % Calculate expected pattern frequency
        Q_i = Q(i, :);
        if any(Q_i == 1)
            n2 = sum(Q_i == 2);
            xi_i = prod(1 - xi(Q_i == 2)) * prod(xi(Q_i ~= 2));
            if n2 == 0
                y_i = x_T(bi2de(Q_i, 'left-msb')) * xi_i;
            else
                U_i = de2bi(0:(2^n2 - 1), n2, 'left-msb');
                y_ij = zeros(2^n2, 1);
                for j = 1:2^n2
                    P_ij = Q_i;
                    P_ij(P_ij == 2) = U_i(j, :);
                    if ~all(P_ij == 0)
                        y_ij(j) = x_T(bi2de(P_ij, 'left-msb'));
                    end
                end
                y_i = sum(y_ij) * xi_i;
            end
            intLogLikExp_i(i) = n_T(i) * log(y_i);
            logLikExp_i(i) = n_T(i) * log(state.lambda * y_i);
            y_T(i) = y_i;
        end
    end
    assertEqual(testCase, sum(x_T) - C, sum(y_T), 'AbsTol', 1e-10);
    intLogLikExp = sum(intLogLikExp_i) - sum(n_T) * log(sum(y_T));
    logLikExp = sum(logLikExp_i) - state.lambda * sum(y_T);
    assertEqual(testCase, intLogLikObs, intLogLikExp, 'AbsTol', 1e-10);
    assertEqual(testCase, logLikObs, logLikExp, 'AbsTol', 1e-10);
end




% function Data3Missing0Registration1(testCase)
% % Test specific code
% end
%
% function Data3Missing1Registration0(testCase)
% % Test specific code
% end
%
% function Data3Missing1Registration1(testCase)
% % Test specific code
% end

% Optional file fixtures
function setupOnce(testCase)
    % Global variables
    GlobalSwitches;
    GlobalValues;
    % Shuffle rng away from start-up seed
    rng('shuffle');
    % Clear persistent variables
    clear patternCounts;
end

function teardownOnce(testCase)
    % Clear persistent variables
    clear logLkd2 patternCounts;
end

% %% Optional fresh fixtures
% function setup(testCase)
% % open a figure, for example
% end
%
% function teardown(testCase)
% % close figure, for example
% end

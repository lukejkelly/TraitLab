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
    clear patternCounts;
    [intLogLikObs, logLikObs] = patternCounts(state, rl, x_T, borPars);
    intLogLikExp = sum(n_T .* log(x_T)) - sum(n_T) * log(sum(x_T) - C);
    logLikExp = sum(n_T .* log(x_T * state.lambda)) ...
                - state.lambda * (sum(x_T) - C);
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

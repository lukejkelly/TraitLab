function tests = discreteUniformLogProbTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    n_i = 1e3;
    n_j = 1e2;
    lpObs = nan(n_i, n_j);
    lpExp = -inf(n_i, n_j);
    for i = 1:n_i
        n = 1 + floor(rand * n_j);
        s = randsample(n_j, n);
        lpObs(i, :) = arrayfun(@(x) discreteUniformLogProb(x, s), 1:n_j);
        lpExp(i, s) = -log(n);
    end
    assertEqual(testCase, lpObs, lpExp);
end

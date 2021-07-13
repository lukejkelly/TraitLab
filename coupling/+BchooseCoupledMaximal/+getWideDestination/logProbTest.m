function tests = logProbTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    n_i = 20;
    for i = 1:n_i
        r = 1:i;
        u = 1 + floor(rand * i);
        v = randsample(i, u);

        lObs = arrayfun(@(x) BchooseCoupledMaximal.getWideDestination.logProb(x, v), r);

        lExp = -inf(size(lObs));
        lExp(v) = -log(u);

        assertEqual(testCase, lObs, lExp);
        assertEqual(testCase, sum(exp(lObs)), 1, 'AbsTol', 1e-12);
    end
end

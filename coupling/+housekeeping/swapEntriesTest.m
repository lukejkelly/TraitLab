function tests = swapEntriesTest
    tests = functiontests(localfunctions);
end

function simpleTest(testCase)
    n = 20;
    x = rand(1, n);
    for i = 1:10
        j = randsample(n, 2);
        xObs = housekeeping.swapEntries(x, j(1), j(2));
        j_min = min(j);
        j_max = max(j);
        xExp = x([1:(j_min - 1), j_max, (j_min + 1):(j_max - 1), j_min, ...
                  (j_max + 1):end]);
        assertEqual(testCase, xObs, xExp);
    end
end

function tests = getOverallCountsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    cat_x = struct('j', 1, 'h', 3, 'p', 0);
    cat_y = struct('j', 2, 'h', 0, 'p', 0);
    catObs = BcatsCoupled.getOverallCounts(cat_x, cat_y);
    catExp = struct('j', 2, 'h', 3, 'p', 0);
    assertEqual(testCase, catObs, catExp);
end

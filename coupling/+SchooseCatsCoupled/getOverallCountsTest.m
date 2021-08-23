function tests = getOverallCountsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    cat_x = struct('i', 1, 'j', 3, 'k', 0);
    cat_y = struct('i', 2, 'j', 0, 'k', 0);
    catObs = SchooseCatsCoupled.getOverallCounts(cat_x, cat_y);
    catExp = struct('i', 2, 'j', 3, 'k', 0);
    assertEqual(testCase, catObs, catExp);
end

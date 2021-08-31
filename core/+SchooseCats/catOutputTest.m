function tests = catOutputTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    catn = ceil(rand(1, 3) * 10);
    catObs = SchooseCats.catOutput(catn);
    assertEqual(testCase, catObs.i, catn(1));
    assertEqual(testCase, catObs.j, catn(2));
    assertEqual(testCase, catObs.k, catn(3));
end

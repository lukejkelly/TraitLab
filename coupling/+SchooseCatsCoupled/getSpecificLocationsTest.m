function tests = getSpecificLocationsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    cat_x = struct('i', poissrnd(3), 'j', poissrnd(4), 'k', poissrnd(5));
    cat_y = struct('i', poissrnd(4), 'j', poissrnd(5), 'k', poissrnd(3));
    cat = SchooseCatsCoupled.getOverallCounts(cat_x, cat_y);
    loc = Bcats.newCatLocations(cat);
    loc_xObs = SchooseCatsCoupled.getSpecificLocations(cat_x, loc);
    loc_yObs = SchooseCatsCoupled.getSpecificLocations(cat_y, loc);

    assertEqual(testCase, loc_xObs.i, loc.i(1:cat_x.i));
    assertEqual(testCase, loc_xObs.j, loc.j(1:cat_x.j));
    assertEqual(testCase, loc_xObs.k, loc.k(1:cat_x.k));

    assertEqual(testCase, loc_yObs.i, loc.i(1:cat_y.i));
    assertEqual(testCase, loc_yObs.j, loc.j(1:cat_y.j));
    assertEqual(testCase, loc_yObs.k, loc.k(1:cat_y.k));

    assertEqual(testCase, ...
                loc_xObs.i(1:min(cat_x.i, cat_y.i)), ...
                loc_yObs.i(1:min(cat_x.i, cat_y.i)));
    assertEqual(testCase, ...
                loc_xObs.j(1:min(cat_x.j, cat_y.j)), ...
                loc_yObs.j(1:min(cat_x.j, cat_y.j)));
    assertEqual(testCase, ...
                loc_xObs.k(1:min(cat_x.k, cat_y.k)), ...
                loc_yObs.k(1:min(cat_x.k, cat_y.k)));

end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global BORROWING
    BORROWING = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function tests = getSpecificLocationsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    cat_x = struct('j', poissrnd(3), 'h', poissrnd(3), 'p', poissrnd(3));
    cat_y = struct('j', poissrnd(3), 'h', poissrnd(3), 'p', poissrnd(3));
    cat = BcatsCoupled.getOverallCounts(cat_x, cat_y);
    loc = Bcats.newCatLocations(cat);
    loc_xObs = BcatsCoupled.getSpecificLocations(cat_x, loc);
    loc_yObs = BcatsCoupled.getSpecificLocations(cat_y, loc);

    assertEqual(testCase, loc_xObs.j, loc.j(1:cat_x.j));
    assertEqual(testCase, loc_xObs.h, loc.h(1:cat_x.h));
    assertEqual(testCase, loc_xObs.p, loc.p(1:cat_x.p));

    assertEqual(testCase, loc_yObs.j, loc.j(1:cat_y.j));
    assertEqual(testCase, loc_yObs.h, loc.h(1:cat_y.h));
    assertEqual(testCase, loc_yObs.p, loc.p(1:cat_y.p));

    assertEqual(testCase, ...
                loc_xObs.j(1:min(cat_x.j, cat_y.j)), ...
                loc_yObs.j(1:min(cat_x.j, cat_y.j)));
    assertEqual(testCase, ...
                loc_xObs.h(1:min(cat_x.h, cat_y.h)), ...
                loc_yObs.h(1:min(cat_x.h, cat_y.h)));
    assertEqual(testCase, ...
                loc_xObs.p(1:min(cat_x.p, cat_y.p)), ...
                loc_yObs.p(1:min(cat_x.p, cat_y.p)));

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

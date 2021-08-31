function tests = newCatLocationsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    global BORROWING
    cat = struct();
    for i = 0:10
        cat.(sprintf('a%i', i)) = i;
    end

    BORROWING = 0;
    loc0 = Bcats.newCatLocations(cat);
    assertEmpty(testCase, loc0);

    BORROWING = 1;
    loc1 = Bcats.newCatLocations(cat);

    assertEmpty(testCase, loc1.a0);
    for i = 1:10
        assertSize(testCase, loc1.(sprintf('a%i', i)), [1, i]);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

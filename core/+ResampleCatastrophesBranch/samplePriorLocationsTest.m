function tests = samplePriorLocationsTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % Simulate catastrophes from prior, compare with relevant part of LogPrior
    global BORROWING
    for BORROWING = 0:1
        for cat = 0:10
            for rep = 1:10
                catloc = ResampleCatastrophesBranch.samplePriorLocations(cat);

                if BORROWING
                    assertEqual(testCase, size(catloc), [cat > 0, cat]);
                    if cat > 0
                        assertGreaterThan(testCase, catloc, 0);
                        assertLessThan(testCase, catloc, 1);
                    end
                else
                    assertEmpty(testCase, catloc);
                end
            end
        end
    end
end

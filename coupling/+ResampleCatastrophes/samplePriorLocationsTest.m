function tests = samplePriorLocationsTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % Simulate catastrophes from prior, compare with relevant part of LogPrior
    global BORROWING
    for BORROWING = 0:1
        for i = 5:15
            cat = round(5 * rand(i, 1));
            catloc = ResampleCatastrophes.samplePriorLocations(cat);

            if BORROWING
                assertEqual(testCase, size(catloc), size(cat));
                j = cat > 0;
                assertEqual(testCase, cellfun(@(x) size(x, 1), catloc(j)), ...
                            ones(sum(j), 1));
                assertEqual(testCase, cellfun(@(x) size(x, 1), catloc(~j)), ...
                            zeros(sum(~j), 1));
                assertEqual(testCase, cellfun(@(x) size(x, 2), catloc), cat);
            else
                assertEmpty(testCase, catloc);
            end
        end
    end
end

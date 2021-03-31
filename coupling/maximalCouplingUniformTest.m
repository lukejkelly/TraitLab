function tests = maximalCouplingUniformTest
    disp('Only simple tests, distribution checked by other coupling unit tests');
    tests = functiontests(localfunctions);
end

function fullOverlapTest(testCase)
    n = 1e3;
    x = zeros(n, 1);
    y = zeros(n, 1);
    for i = 1:n
        [x(i), y(i)] = maximalCouplingUniform(0, 1, 0, 1);
    end
    assertEqual(testCase, x, y)
end

function noOverlapTest(testCase)
    n = 1e3;
    x = zeros(n, 1);
    y = zeros(n, 1);
    for i = 1:n
        [x(i), y(i)] = maximalCouplingUniform(0, 1, 1, 2);
    end
    assertNotEqual(testCase, x, y)
end

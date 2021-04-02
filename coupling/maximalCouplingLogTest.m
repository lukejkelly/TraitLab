function tests = maximalCouplingLogTest
    disp('Only simple tests, distribution checked by other coupling unit tests');
    tests = functiontests(localfunctions);
end

function fullOverlapTest(testCase)
    n = 1e3;
    x = zeros(n, 1);
    y = zeros(n, 1);
    for i = 1:n
        [x(i), y(i)] = maximalCouplingLog(...
            @() rand(), @(x) log(double(0 <= x && x <= 1)), ...
            @() rand(), @(y) log(double(0 <= y && y <= 1)));
    end
    assertEqual(testCase, x, y)
end

function noOverlapTest(testCase)
    n = 1e3;
    x = zeros(n, 1);
    y = zeros(n, 1);
    for i = 1:n
        [x(i), y(i)] = maximalCouplingLog(...
            @() rand(), @(x) log(double(0 <= x && x <= 1)), ...
            @() 1 + rand(), @(y) log(double(1 <= y && y <= 2)));
    end
    assertNotEqual(testCase, x, y)
end

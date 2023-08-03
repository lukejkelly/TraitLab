function tests = lpdfGTest
    tests = functiontests(localfunctions);
end

function valueTest(testCase)
    % compute log PDF of each x[i] ~ G(a[i], b[i]), compare to R dgamma
    a = [4.966, 3.924, 2.692, 1.394, 1.407, 4.514, 4.360, 1.452, 0.478, 3.163];
    b = [0.696, 3.689, 2.541, 4.821, 1.469, 1.177, 0.522, 2.200, 2.083, 4.929];
    x = [5.613, 1.073, 1.776, 0.447, 0.357, 2.595, 5.693, 0.619, 0.404, 0.835];

    lObs = arrayfun(@(i) lpdfG(x(i), a(i), b(i)), 1:length(x));
    lExp = [...
        -1.9916399, -0.3271819, -1.4590080, -0.1602723, -0.2825267, ...
        -1.4409758, -2.2240516, -0.3123204, -0.6344488, -0.3090389 ...
    ];
    assertEqual(testCase, lObs, lExp, 'AbsTol', 1e-6);
end

function tests = validTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Compare new and old implementations
    for h = 1:1e2
        L = 8 + ceil(rand * 8);
        N = 2 * L - 1;

        i = ceil(N * rand);
        n = 5 + ceil(rand * (N - 5));
        r = randsample(N, n)';

        state = unitTests.dummyState(ExpTree(L, 1e-4));
        s = state.tree;

        vObs = BchooseCoupled.getWideDestination.valid(i, r, s);

        uExp = ones(1, n);
        iT = s(i).time;
        for j = r
            k = s(j).parent;
            if (s(k).time <= iT || i == j || i == k)
                uExp(ismember(r, j)) = 0;
            end
        end
        vExp = r(uExp == 1);

        assertEqual(testCase, sort(vObs), sort(vExp));
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

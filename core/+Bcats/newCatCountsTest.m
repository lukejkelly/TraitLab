function tests = newCatCountsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % compare independent(-ish) implementations
    global ROOT
    M = 20;
    for m = 1:M
        L = 8 + ceil(rand * 8);
        state = dummyState(L);
        s = state.tree;
        r = state.root;
        v = find([s.type] < ROOT);

        for i = v
            for j = setdiff(v, i)
                if j == r
                    jn = poissrnd(1);
                else
                    jn = binornd(state.cat(j), 2 / 3);
                end
                [ncatObs, catObs] = Bcats.newCatCounts(state, i, j, jn);

                [p, ~, h] = Bcats.pqh(state, i);
                catExp = struct('j', jn, 'h', [], 'p', []);
                switch r
                case p
                    ncatExp = state.ncat - state.cat(h);
                    catExp.h = 0;
                    catExp.p = state.cat(j) - jn;
                case j
                    ncatExp = state.ncat + jn;
                    catExp.h = sum(state.cat([h, p]));
                    catExp.p = 0;
                otherwise
                    ncatExp = state.ncat;
                    catExp.h = sum(state.cat([h, p]));
                    catExp.p = state.cat(j) - jn;
                end

                assertEqual(testCase, ncatObs, ncatExp);
                assertEqual(testCase, catObs, catExp);
            end
        end
    end
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, rand * 1e-2));
    for c = 1:poissrnd(3)
        state = AddCat(state);
    end
end

% Setup and teardown functions
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
    global MCMCCAT BORROWING;
    MCMCCAT = 1;
    BORROWING = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

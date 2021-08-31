function tests = updateCatsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % compare independent(-ish) implementations
    global ROOT BORROWING
    L = 8;
    for BORROWING = [0, 1]
        state = dummyState(L);
        s = state.tree;
        r = state.root;
        v = find([s.type] < ROOT);
        for i = v
            [p, ~, h] = Bcats.pqh(state, i);
            for j = setdiff(v, [i, p, h])
                if j == r
                    jn = poissrnd(1);
                else
                    jn = binornd(state.cat(j), 2 / 3);
                end
                [ncat, cat] = Bcats.newCatCounts(state, i, j, jn);
                loc = Bcats.newCatLocations(cat);
                nstate = Bcats.updateCats(state, i, j, ncat, cat, loc);

                assertEqual(testCase, nstate.ncat, ncat);
                if p == r
                    assertLessThanOrEqual(testCase, nstate.ncat, state.ncat);
                elseif j == r
                    assertLessThanOrEqual(testCase, nstate.ncat, state.ncat);
                else
                    assertEqual(testCase, nstate.ncat, state.ncat);
                end

                assertEqual(testCase, nstate.cat(j), cat.j);
                assertEqual(testCase, nstate.cat(h), cat.h);
                assertEqual(testCase, nstate.cat(p), cat.p);

                if BORROWING
                    assertEqual(testCase, nstate.tree(j).catloc, sort(loc.j));
                    assertEqual(testCase, nstate.tree(h).catloc, sort(loc.h));
                    assertEqual(testCase, nstate.tree(p).catloc, sort(loc.p));
                else
                    assertEmpty(testCase, nstate.tree(j).catloc);
                    assertEmpty(testCase, nstate.tree(h).catloc);
                    assertEmpty(testCase, nstate.tree(p).catloc);
                end
                for k = setdiff(v, [j, h, p])
                    assertEqual(testCase, nstate.cat(k), state.cat(k));
                    assertEqual(testCase, nstate.tree(k).catloc, state.tree(k).catloc);
                end
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
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

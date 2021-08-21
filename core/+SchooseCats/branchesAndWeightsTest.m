function tests = branchesAndWeightsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    for L = 6:10
        state = dummyState(L);
        r = state.root;
        s = state.tree;
        for j = 1:length(state.nodes)
            i = state.nodes(j);
            [iT, kT, jT, a, b] = SchooseTime.nodeTimesAndRanges(i, state);
            if i == r
                newage = a + rand * (b - a);
            else
                newage = jT + rand * (kT - jT);
            end

            [catcObs, ncatObs, wcObs, wnObs] = SchooseCats.branchesAndWeights(...
                state, i, newage);

            bExp = [i, s(i).child];
            catcExp = state.cat(bExp);
            assertEqual(testCase, catcObs, catcExp);

            ncatExp = sum(catcExp);
            assertEqual(testCase, ncatObs, ncatExp);

            pT = s(s(i).parent).time;
            cT = [s(s(i).child).time];
            if i == r
                vcExp = [0, iT  - cT];
            else
                vcExp = [pT - iT, iT - cT];
            end
            wcExp = vcExp ./ sum(vcExp);
            assertEqual(testCase, wcObs, wcExp, 'AbsTol', 1e-10);

            if i == r
                vnExp = vcExp + (iT - newage) * [0, -1, -1];
            else
                vnExp = vcExp + (iT - newage) * [+1, -1, -1];
            end
            wnExp = vnExp ./ sum(vnExp);
            assertEqual(testCase, wnObs, wnExp, 'AbsTol', 1e-10);
        end
    end
end

function setupOnce(testCase)
    global MCMCCAT
    unitTests.setupOnce(testCase);
    MCMCCAT = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
    close;
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for i = 1:(3 + poissrnd(3))
        state = AddCat(state);
    end
end

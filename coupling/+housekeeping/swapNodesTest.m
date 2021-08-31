function tests = swapNodesTest
    tests = functiontests(localfunctions);
end

function orderTest(testCase)
    for L = randsample(2:15, 1e2, true)
        s = ExpTree(L, 1);
        jk = randsample(2 * L, 2);
        sObs1 = housekeeping.swapNodes(s, jk(1), jk(2));
        sObs2 = housekeeping.swapNodes(s, jk(2), jk(1));
        assertEqual(testCase, sObs1, sObs2);
    end
end

function swapTest(testCase)
    % Apart from sibling information when j and k are direct relatives
    for L = randsample(2:15, 1e2, true)
        s = ExpTree(L, 1);
        jk = randsample(2 * L, 2);
        j = jk(1);
        k = jk(2);
        sObs = housekeeping.swapNodes(s, j, k);
        sExp = swap(s, j, k);
        if any([s(jk).parent] == [k, j])
            if s(k).parent == j
                [j, k] = deal(k, j);
            end
            if s(j).sibling == 1
                assertEqual(testCase, sObs, sExp);
            else
                for fn = {'Name', 'parent', 'type', 'time'}
                    assertEqual(testCase, {sObs.(fn{:})}, {sExp.(fn{:})});
                end
                assertEqual(testCase, sObs(j).child, fliplr(sExp(j).child));
                assertEqual(testCase, sObs(k).sibling, 2);
                assertEqual(testCase, sExp(k).sibling, 1);
                c = s(k).child(s(k).child ~= j);
                if ~isempty(c)
                    assertEqual(testCase, sObs(c).sibling, 1);
                    assertEqual(testCase, sExp(c).sibling, 2);
                end
                assertEqual(testCase, equaltrees(sObs, sExp), 1);
            end
        else
            assertEqual(testCase, sObs, sExp);
        end
    end
end

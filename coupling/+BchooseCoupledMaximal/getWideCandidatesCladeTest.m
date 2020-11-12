function tests = getWideCandidatesCladeTest
    % Unit-testing various maximal coupling functions
    tests = functiontests(localfunctions);
end

function cladesYesTest(testCase)
    % Create struct of elegible destinations for SPR moves <pa(i), i)> grafted
    % onto <pa(a), a> so pa(a) becomes a child of pa(a) and sibling of a
    % For each node i we want the nodes a such that
    %     cl(pa(a)) = cl[pa(i)] SUBSET cl(a)
    % In other words, pa(i) is as constrained as pa(a) and least as constrained
    % as node a, which also cannot be the root
    % We ignore time information for now
    global LEAF ANST

    s = BchooseCoupledMaximal.state10('Yes');

    rExp = cell(size(s));
    rExp{1} = [];
    rExp{2} = [2, 3, 4, 13, 14, 15, 16, 17, 18, 19];
    rExp{3} = rExp{2};
    rExp{4} = rExp{2};
    rExp{5} = [5, 6, 7, 8, 9, 10];
    rExp{6} = rExp{5};
    rExp{7} = rExp{5};
    rExp{8} = rExp{5};
    rExp{9} = rExp{5};
    rExp{10} = rExp{5};
    rExp{11} = [11, 12];
    rExp{12} = [11, 12];
    rExp{13} = rExp{2};
    rExp{14} = rExp{2};
    rExp{15} = rExp{2};
    rExp{16} = rExp{2};
    rExp{17} = rExp{2};
    rExp{18} = rExp{2};
    rExp{19} = rExp{2};
    rExp{20} = [];

    for i = find(ismember([s.type], [LEAF, ANST]))
        rObs_i = BchooseCoupledMaximal.getWideCandidatesClade(i, s);
        assertEqual(testCase, rObs_i, rExp{i});
    end
    warning('Unlike mcmc with no clade prior, we cannot SPR to <Adam, root>??');
end

function cladesNoTest(testCase)
    % We can regraft any branch onto any other one, indexed by the child in
    % each case
    global LEAF ANST

    s = BchooseCoupledMaximal.state10('No');

    rExp = find(ismember([s.type], [LEAF, ANST]));
    for i = find(ismember([s.type], [LEAF, ANST]))
        rObs_i = BchooseCoupledMaximal.getWideCandidatesClade(i, s);
        assertEqual(testCase, rObs_i, rExp);
    end
    warning('Unlike mcmc with no clade prior, we cannot SPR to <Adam, root>??');
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

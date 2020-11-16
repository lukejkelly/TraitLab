function tests = getWideDestinationTest
    tests = functiontests(localfunctions);
end

function cladesYesTest(testCase)
    global ANST

    s = BchooseCoupledMaximal.state10('Yes');
    for i = find([s.type] <= ANST)
        r_i = BchooseCoupledMaximal.getWideCandidatesClade(i, s);
        N_i = length(r_i);
        for iPerm = 1:N_i
            rPerm = r_i(randperm(N_i));
            for iShift = 1:N_i
                rShift = circshift(rPerm, iShift);

                [jObs, kObs, FAILObs] ...
                    = BchooseCoupledMaximal.getWideDestination(i, rShift, ...
                                                               N_i, s);
                if N_i > 4
                    jExp = rShift(find(arrayfun(...
                        @(a) s(s(a).parent).time > s(i).time ...
                             && ~ismember(i, [a, s(a).parent]), ...
                        rShift), 1));
                    kExp = s(jExp).parent;
                    FAILExp = 0;
                    if ismember(s(i).parent, [jExp, kExp])
                        FAILExp = 1;
                    end
                else
                    [jExp, kExp] = deal(-1);
                    FAILExp = 1;
                end
                assertEqual(testCase, jObs, jExp);
                assertEqual(testCase, kObs, kExp);
                assertEqual(testCase, FAILObs, FAILExp);
            end
        end
    end

    warning(['The logic here could be greatly simplified and sampling ', ...
             'efficiency improved as we are proposing destinations only ', ...
             'to discard them and thereby reject the move']);
end

function cladesNoTest(testCase)
    global ANST ROOT

    s = BchooseCoupledMaximal.state10('No');
    r = find([s.type] <= ROOT);
    N = length(r);
    for i = find([s.type] <= ANST)
        for iPerm = 1:N
            rPerm = r(randperm(N));
            for iShift = 1:N
                rShift = circshift(rPerm, iShift);

                [jObs, kObs, FAILObs] ...
                    = BchooseCoupledMaximal.getWideDestination(i, rShift, ...
                                                               N, s);
                if N > 4
                    jExp = rShift(find(arrayfun(...
                        @(a) s(s(a).parent).time > s(i).time ...
                             && ~ismember(i, [a, s(a).parent]), ...
                        rShift), 1));
                    kExp = s(jExp).parent;
                    FAILExp = 0;
                    if ismember(s(i).parent, [jExp, kExp])
                        FAILExp = 1;
                    end
                else
                    [jExp, kExp] = deal(-1);
                    FAILExp = 1;
                end
                assertEqual(testCase, jObs, jExp);
                assertEqual(testCase, kObs, kExp);
                assertEqual(testCase, FAILObs, FAILExp);
            end
        end
    end

    warning(['The logic here could be greatly simplified and sampling ', ...
             'efficiency improved as we are proposing destinations only ', ...
             'to discard them and thereby reject the move']);
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

function tests = getUniformParametersTest
    tests = functiontests(localfunctions);
end

function iPRootYesTest(testCase)
    % Moving the root
    global ROOT OTHER
    THETA = testCase.TestData.THETA;

    s = BchooseCoupledMaximal.state10a('No');

    for i = find(arrayfun(@(a) ~isempty(s(a).parent) ...
                          && s(s(a).parent).type == ROOT, 1:length(s)))

        iT = s(i).time;
        iP = s(i).parent;
        iPT = s(iP).time;
        CiPT = s(s(iP).child(OTHER(s(i).sibling))).time;
        logq_reverse = log(THETA) - THETA * (iPT - CiPT);

        for j = find(arrayfun(@(a) ~isempty(s(a).parent) ...
                              && i ~= a && i ~= s(a).parent ...
                              && iP ~= a && iP ~= s(a).parent ...
                              && s(s(a).parent).time > iT, 1:length(s)))

            k = s(j).parent;

            [new_minageObs, kTObs, logqObs] ...
                = BchooseCoupledMaximal.sampleCoupling.getUniformParameters(...
                    i, j, k, s, THETA);

            new_minageExp = max(iT, s(j).time);
            kTExp = s(k).time;
            logqExp = logq_reverse + log(kTExp - new_minageExp);

            assertEqual(testCase, new_minageObs, new_minageExp);
            assertEqual(testCase, kTObs, kTExp);
            assertEqual(testCase, logqObs, logqExp, 'RelTol', 1e-12);
        end
    end
end

function iPRootNoTest(testCase)
    % Neither source nor destination includes root
    global ROOT OTHER
    THETA = testCase.TestData.THETA;

    s = BchooseCoupledMaximal.state10a('No');

    for i = find(arrayfun(@(a) ~isempty(s(a).parent) ...
                          && s(s(a).parent).type < ROOT, 1:length(s)))

        iT = s(i).time;
        iP = s(i).parent;
        logq_reverse = -log(s(s(iP).parent).time - max(s(s(iP).child).time));

        for j = find(arrayfun(@(a) ~isempty(s(a).parent) ...
                              && i ~= a && i ~= s(a).parent ...
                              && iP ~= a && iP ~= s(a).parent ...
                              && s(s(a).parent).time > iT, 1:length(s)))

            k = s(j).parent;

            [new_minageObs, kTObs, logqObs] ...
                = BchooseCoupledMaximal.sampleCoupling.getUniformParameters(...
                    i, j, k, s, THETA);

            new_minageExp = max(iT, s(j).time);
            kTExp = s(k).time;
            logqExp = logq_reverse + log(kTExp - new_minageExp);

            assertEqual(testCase, new_minageObs, new_minageExp);
            assertEqual(testCase, kTObs, kTExp);
            assertEqual(testCase, logqObs, logqExp, 'RelTol', 1e-12);
        end
    end
end

function setupOnce(testCase)
    GlobalSwitches;
    GlobalValues;
    testCase.TestData.THETA = 0.001;
end

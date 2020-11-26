function tests = getShiftedExponentialParameterTest
    % Destination branch is always <Adam, Root>
    tests = functiontests(localfunctions);
end

function jRootYesTest(testCase)
    global ROOT
    THETA = 0.01;

    s = BchooseCoupledMaximal.state10('No');

    j = find([s.type] == ROOT);
    jT = s(j).time;
    k = s(j).parent;

    for i = find(arrayfun(@(a) a ~= j && ~isempty(s(a).parent) ...
                          && s(a).parent ~= j, 1:length(s)))
         [newage, logqObs] = BchooseCoupledMaximal.sampleMarginal(i, j, k, ...
                                                                  s, THETA);

         iPPT = s(s(s(i).parent).parent).time;
         iPeCT = max(s(s(s(i).parent).child).time);
         logqExp = THETA * (newage - jT) - log(THETA) - log(iPPT - iPeCT);

         assertEqual(testCase, logqObs, logqExp, 'RelTol', 1e-12);
    end
end

function setupOnce(~)
    GlobalSwitches;
    GlobalValues;
end

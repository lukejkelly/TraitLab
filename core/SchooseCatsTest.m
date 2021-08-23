function tests = SchooseCatsTest
    tests = functiontests(localfunctions);
end

function outputTest(testCase)
    % Only checking outputs as distributions checked elsewhere
    global BORROWING
    for BORROWING = 0:1
        for L = 8:12
            state = dummyState(L);
            r = state.root;
            for j = 1:length(state.nodes)
                i = state.nodes(j);
                [~, kT, jT, a, b] = SchooseTime.nodeTimesAndRanges(i, state);
                if i == r
                    newage = a + rand * (b - a);
                else
                    newage = jT + rand * (kT - jT);
                end
                [cat, loc, ~] = SchooseCats(state, i, newage);
                iC = [state.tree(i).child];
                assertEqual(testCase, ...
                            cat.i + cat.j + cat.k, ...
                            sum(state.cat([i, iC])));
                if BORROWING
                    assertEqual(testCase, ...
                                structfun(@length, loc), ...
                                [cat.i; cat.j; cat.k]);
                else
                    assertEmpty(testCase, loc);
                end
            end
        end
    end
end

function manualTest(testCase)
    % Simple distribution checks
    N = 1e4;
    global BORROWING
    for BORROWING = 0:1
        fprintf('BORROWING = %i\n', BORROWING);

        state = unitTests.dummyState(rnextree('((1:1, 2:1):1, 3:2);'));
        %      ---1---   Root time = 2
        %     |       |
        %   --2--     |  Node 2 time = 1
        %  |     |    |
        % 3(1) 4(2)  5(3)

        state.cat(2:3) = 1;
        if BORROWING
            state.tree(2).catloc = rand;
            state.tree(3).catloc = rand;
        end
        state.ncat = sum(state.cat);

        p1 = (0:2) / 3;
        n1 = [0, 1, 0; 0, 0, 1];
        q1 = mnpdf(n1, p1);
        l1 = log(q1(1)) - log(q1);

        cat1 = zeros(size(n1, 1), 1);
        [loc1, logq1] = deal(nan(N, 1));
        for n = 1:N
            [cat, loc, logq] = SchooseCats(state, 1, state.tree(1).time);
            cat_n = arrayfun(@(s) cat.(s), 'ijk');
            for j = 1:size(n1, 1)
                if isequal(cat_n, n1(j, :))
                    cat1(j) = cat1(j) + 1;
                    if BORROWING
                        loc1(n) = isequal(...
                            cat_n, ...
                            ~cellfun(@isempty, ...
                                     arrayfun(@(s) loc.(s), 'ijk', ...
                                              'UniformOutput', false)));
                    else
                        loc1(n) = isempty(loc);
                    end
                    logq1(n) = logq == l1(j);
                end
            end
        end
        fprintf('i = 1\n');
        fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cat1 / N));
        fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', q1));
        fprintf('Difference: %s\n\n', sprintf('%+-3.2e\t', cat1 / N - q1));
        assertTrue(testCase, all(loc1));
        assertTrue(testCase, all(logq1));

        p2 = ones(1, 3) / 3;
        n2 = [2 * eye(3); [1, 1, 0]; [1, 0, 1]; [0, 1, 1]];
        q2 = mnpdf(n2, p2);
        l2 = log(q2(4)) - log(q2);

        cat2 = zeros(size(n2, 1), 1);
        [loc2, logq2] = deal(nan(N, 1));
        for n = 1:N
            [cat, loc, logq] = SchooseCats(state, 2, state.tree(2).time);
            cat_n = arrayfun(@(s) cat.(s), 'ijk');
            for j = 1:size(n2, 1)
                if isequal(cat_n, n2(j, :))
                    cat2(j) = cat2(j) + 1;
                    if BORROWING
                        loc2(n) = isequal(...
                            cat_n, ...
                            cellfun(@length, ...
                                     arrayfun(@(s) loc.(s), 'ijk', ...
                                              'UniformOutput', false)));
                    else
                        loc2(n) = isempty(loc);
                    end
                    logq2(n) = logq == l2(j);
                end
            end
        end
        fprintf('i = 2\n');
        fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cat2 / N));
        fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', q2));
        fprintf('Difference: %s\n\n', sprintf('%+-3.2e\t', cat2 / N - q2));
        assertTrue(testCase, all(loc2));
        assertTrue(testCase, all(logq2));
    end
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end


function setupOnce(testCase)
    global MCMCCAT
    unitTests.setupOnce(testCase);
    MCMCCAT = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function state = dummyState(L)
    state = unitTests.dummyState(ExpTree(L, 1e-2));
    for i = 1:(5 + poissrnd(5))
        state = AddCat(state);
    end
end

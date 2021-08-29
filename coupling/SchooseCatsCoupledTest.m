function tests = SchooseCatsCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    global BORROWING
    M = 10;
    lambda_max = 5;
    for BORROWING = [0, 1]
        for L = 8:12
            for m = 1:M
                [state, ~] = coupledStates(L);
                for lambda = 0:lambda_max
                    [state_x, state_y] = deal(state);
                    for c = 1:poissrnd(lambda)
                        [state_x, state_y] = AddCatCoupled(state_x, state_y);
                    end
                    [i, newage_x, newage_y, ~, ~] = SchooseTimeCoupled(...
                        state_x, state_y);
                    [cat_x, cat_y, loc_x, loc_y, logq_x_c, logq_y_c] ...
                        = SchooseCatsCoupled(state_x, state_y, i, ...
                                             newage_x, newage_y);
                    assertEqual(testCase, cat_x, cat_y);
                    assertEqual(testCase, loc_x, loc_y);
                    assertEqual(testCase, logq_x_c, logq_y_c);
                end
            end
        end
    end
end

function couplingTest(testCase)
    % Only consider case with same number of catastrophes at each triplet
    global BORROWING
    N = 1e4;
    Nis = 1e5;
    L = 8;
    fprintf('Proportion of matching samples at each internal node\n');
    fprintf('Observed is simple MC average from %.1e trials\n', N);
    fprintf('Expected is IS estimate from %.1e trials\n\n', Nis);

    for BORROWING = [0, 1]
        [cObs, cExp, lObs, lExp] = deal(zeros(1, L - 1));
        [state_x, state_y] = coupledStates(L);
        r = state_x.root;
        for c = 1:poissrnd(5)
            [state_x, state_y] = AddCatCoupled(state_x, state_y);
        end

        for j = 1:length(state_x.nodes)
            iFound = false;
            while ~iFound
                [i, newage_x, ~, cat_x, loc_x] = Schoose(state_x);
                iFound = i == state_x.nodes(j);
            end
            nstate_x = Supdate(state_x, i, newage_x, cat_x, loc_x);

            s_x = nstate_x.tree;
            s_y = state_y.tree;

            % Same indices in both states
            pa = s_x(i).parent;
            ch = s_x(i).child;

            d_x = [s_x(pa).time - newage_x, newage_x - [s_x(ch).time]];

            newage_y = s_y(i).time;
            d_y = [s_y(pa).time - newage_y, newage_y - [s_y(ch).time]];

            if i == r
                d_x(1) = 0;
                d_y(1) = 0;
            end
            f_x = d_x / sum(d_x);
            f_y = d_y / sum(d_y);

            for n = 1:N
                [cat_x, cat_y, loc_x, loc_y, ~, ~] = SchooseCatsCoupled(...
                    nstate_x, state_y, i, newage_x, newage_y);

                if isequal(cat_x, cat_y)
                    cObs(j) = cObs(j) + 1;
                end

                if BORROWING
                    lfun = @(s) isequal(...
                        length(intersect(loc_x.(s), loc_y.(s))), ...
                        min(cat_x.(s), cat_y.(s)));
                    if all(arrayfun(@(s) lfun(s), 'ijk'))
                        lObs(j) = lObs(j) + 1;
                    end
                else
                    if isempty(loc_x) && isempty(loc_y)
                        lObs(j) = lObs(j) + 1;
                    end
                end
            end

            % Importance sampling estimate of overlap
            ncat = sum(arrayfun(@(s) cat_x.(s), 'ijk'));
            f_u = ones(1, 3) / 3;
            z = mnrnd(ncat, f_u, Nis);
            cExp(j) = mean(min(mnpdf(z, f_x), mnpdf(z, f_y)) ./ mnpdf(z, f_u));

            % Intersection of locations should always match
            lExp(j) = N;

            if abs(cObs(j) / N - cExp(j)) > 5e-2
                keyboard;
            end
        end
        assertEqual(testCase, lObs, lExp);

        fprintf('BORROWING = %i\n', BORROWING);
        fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs / N));
        fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
        fprintf('Difference: %s\n\n', sprintf('%+-3.2e\t', cObs / N - cExp));
    end
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function housekeptTest(testCase)
    % Pick a node at random, in most cases there probably won't be an overlap
    % in catastrophe counts
    global BORROWING
    M = 8;
    N = 1e4;
    Nis = 1e5;
    fprintf('Proportion of matching samples\n');
    fprintf('Observed is simple MC average from %.1e trials\n', N);
    fprintf('Expected is IS estimate from %.1e trials\n\n', Nis);

    for BORROWING = [0, 1]
        [cObs, cExp, lObs, lExp] = deal(zeros(1, M));
        for m = 1:8
            L = 5 + m;
            [state_x, state_y] = housekeptStates(L);
            r = state_x.root;
            for c = 1:poissrnd(3)
                [state_x, state_y] = AddCatCoupled(state_x, state_y);
            end

            [i, newage_x, ~, cat_x, loc_x] = Schoose(state_x);
            nstate_x = Supdate(state_x, i, newage_x, cat_x, loc_x);

            s_x = nstate_x.tree;
            s_y = state_y.tree;

            % Same indices in both states
            p_x = s_x(i).parent;
            c_x = s_x(i).child;

            d_x = [s_x(p_x).time - newage_x, newage_x - [s_x(c_x).time]];

            p_y = s_y(i).parent;
            c_y = s_y(i).child;
            newage_y = s_y(i).time;
            d_y = [s_y(p_y).time - newage_y, newage_y - [s_y(c_y).time]];

            if i == r
                d_x(1) = 0;
                d_y(1) = 0;
            end
            f_x = d_x / sum(d_x);
            f_y = d_y / sum(d_y);

            for n = 1:N
                [cat_x, cat_y, loc_x, loc_y, ~, ~] = SchooseCatsCoupled(...
                    nstate_x, state_y, i, newage_x, newage_y);

                if isequal(cat_x, cat_y)
                    cObs(m) = cObs(m) + 1;
                end

                if BORROWING
                    lfun = @(s) isequal(...
                        length(intersect(loc_x.(s), loc_y.(s))), ...
                        min(cat_x.(s), cat_y.(s)));
                    if all(arrayfun(@(s) lfun(s), 'ijk'))
                        lObs(m) = lObs(m) + 1;
                    end
                else
                    if isempty(loc_x) && isempty(loc_y)
                        lObs(m) = lObs(m) + 1;
                    end
                end
            end

            % Importance sampling estimate of overlap
            ncat_x = sum(arrayfun(@(s) cat_x.(s), 'ijk'));
            ncat_y = sum(arrayfun(@(s) cat_y.(s), 'ijk'));
            if ncat_x == ncat_y
                f_u = ones(1, 3) / 3;
                z = mnrnd(ncat_x, f_u, Nis);
                cExp(m) = mean(min(mnpdf(z, f_x), mnpdf(z, f_y)) ...
                               ./ mnpdf(z, f_u));
            else
                cExp(m) = 0;
            end

            % Intersection of locations should always match
            lExp(m) = N;

            if abs(cObs(m) / N - cExp(m)) > 5e-2
                keyboard;
            end
        end
        assertEqual(testCase, lObs, lExp);

        fprintf('BORROWING = %i\n', BORROWING);
        fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs / N));
        fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
        fprintf('Difference: %s\n\n', sprintf('%+-3.2e\t', cObs / N - cExp));
    end
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end

function manualTest(testCase)
    % Simple distribution checks
    global BORROWING
    N = 1e4;
    fprintf('Coupling proportion in %i trials at single node\n\n', BORROWING);
    for BORROWING = 0:1
        fprintf('BORROWING = %i\n', BORROWING);
        state_x = unitTests.dummyState(rnextree('((1:1, 2:1):1, 3:2);'));
        %      ---1---   Root time = 2
        %     |       |
        %   --2--     |  Node 2 time = 1
        %  |     |    |
        % 3(1) 4(2)  5(3)

        % Identical to x except node 2 occurs at time 0.5
        state_y = unitTests.dummyState(rnextree('((1:0.5, 2:0.5):1.5, 3:2);'));

        state_x.cat(2:3) = 1;
        state_y.cat(3:4) = 1;
        if BORROWING
            state_x.tree(2).catloc = rand;
            state_x.tree(3).catloc = rand;
            state_y.tree(3).catloc = rand;
            state_y.tree(4).catloc = rand;
        end
        state_x.ncat = sum(state_x.cat);
        state_y.ncat = sum(state_y.cat);

        i = 2;
        newage_x = state_x.tree(i).time;
        newage_y = state_y.tree(i).time;

        p_x = ones(1, 3) / 3;
        p_y = [3, 1, 1] / 5;

        n2 = [2 * eye(3); [1, 1, 0]; [1, 0, 1]; [0, 1, 1]];
        q_x = mnpdf(n2, p_x);
        q_y = mnpdf(n2, p_y);

        cExp = sum(min(q_x, q_y));

        l_x = log(q_x(4)) - log(q_x);
        l_y = log(q_y(6)) - log(q_y);

        cObs = 0;
        [loc, logq] = deal(nan(N, 1));
        for n = 1:N
            [cat_x, cat_y, loc_x, loc_y, logq_x, logq_y] ...
                = SchooseCatsCoupled(state_x, state_y, i, newage_x, newage_y);

            c_x = arrayfun(@(s) cat_x.(s), 'ijk');
            c_y = arrayfun(@(s) cat_y.(s), 'ijk');
            if all(c_x == c_y)
                cObs = cObs + 1;
            end

            if BORROWING
                if all(c_x == c_y)
                    loc(n) = isequal(loc_x, loc_y);
                else
                    loc(n) = ~isequal(loc_x, loc_y);
                end
            else
                loc(n) = isempty(loc_x) && isempty(loc_y);
            end

            j_x = all(n2 == c_x, 2);
            j_y = all(n2 == c_y, 2);
            logq(n) = (logq_x == l_x(j_x)) && (logq_y == l_y(j_y));
        end
        fprintf('Observed:   %s\n', sprintf(' %-3.2e\t', cObs / N));
        fprintf('Expected:   %s\n', sprintf(' %-3.2e\t', cExp));
        fprintf('Difference: %s\n\n', sprintf('%+-3.2e\t', cObs / N - cExp));
        assertTrue(testCase, all(loc));
        assertTrue(testCase, all(logq));
    end
    v = input('Are these proportions the same? Reply 1 for yes... ');
    assertTrue(testCase, v == 1);
end


%%%
function [state_x, state_y] = coupledStates(L)
    [state_x, state_y] = unitTests.coupledStates(L, 5e-3 * rand);
end

function [state_x, state_y] = housekeptStates(L)
    [state_x, state_y] = unitTests.housekeptStates(L, 5e-3 * rand);
end

% Setup and teardown functions
function setupOnce(testCase)
    global MCMCCAT
    unitTests.setupOnce(testCase);
    MCMCCAT = 1;
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

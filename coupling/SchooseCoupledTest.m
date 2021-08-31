function tests = SchooseCoupledTest
    tests = functiontests(localfunctions);
end

function coupledTest(testCase)
    % Only checking outputs as distributions checked elsewhere
    global BORROWING MCMCCAT
    M = 5;
    N = 1e2;
    for BORROWING = 0:1
        for MCMCCAT = 0:1
            [newage, logq, cat, loc] = deal(nan(M, N));
            for m = 1:M
                L = 5 + m;
                for n = 1:N
                    [state_x, state_y] = coupledStates(L);
                    [~, newage_x, newage_y, logq_x, logq_y, cat_x, cat_y, ...
                     loc_x, loc_y] = SchooseCoupled(state_x, state_y);

                     newage(m, n) = isequal(newage_x, newage_y);
                     logq(m, n) = isequal(logq_x, logq_y);
                     if MCMCCAT
                         cat(m, n) = isequal(cat_x, cat_y);
                         if BORROWING
                             loc(m, n) = isequal(loc_x, loc_y);
                         else
                             loc(m, n) = isempty(loc_x) && isempty(loc_y);
                         end
                     else
                         cat(m, n) = isempty(cat_x) && isempty(cat_y);
                         loc(m, n) = isempty(loc_x) && isempty(loc_y);
                     end
                 end
             end
            assertTrue(testCase, all(newage, 'all'));
            assertTrue(testCase, all(logq, 'all'));
            assertTrue(testCase, all(cat, 'all'));
            assertTrue(testCase, all(loc, 'all'));
        end
    end
end

%%%
function setupOnce(testCase)
    unitTests.setupOnce(testCase);
end

function teardownOnce(testCase)
    unitTests.teardownOnce(testCase);
end

function [state_x, state_y] = coupledStates(L)
    global MCMCCAT
    [state_x, state_y] = unitTests.coupledStates(L, 1e-2);
    if MCMCCAT
        for i = 1:(5 + poissrnd(5))
            [state_x, state_y] = AddCatCoupled(state_x, state_y);
        end
    end
end

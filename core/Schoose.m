function [i, newage, logq, cat, loc] = Schoose(state)
    % Sample new node time and possibly update catastrophes
    global MCMCCAT

    % Previous Schoose function
    [i, newage, logq_t] = SchooseTime(state);

    % Shuffle cats around node that moved
    if MCMCCAT
        [cat, loc, logq_c] = SchooseCats(state, i, newage);
    else
        [cat, loc] = deal([]);
        logq_c = 0;
    end

    logq = logq_t + logq_c;
end

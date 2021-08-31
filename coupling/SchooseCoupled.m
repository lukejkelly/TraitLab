function [i, newage_x, newage_y, logq_x, logq_y, cat_x, cat_y, loc_x, loc_y] ...
        = SchooseCoupled(state_x, state_y)
    % Resample node time and shuffle catastrophes on affected branches
    global MCMCCAT

    % Previous SchooseCoupled function
    [i, newage_x, newage_y, logq_x_t, logq_y_t] = SchooseTimeCoupled(...
        state_x, state_y);

    % Shuffle cats around nodes that moved
    if MCMCCAT
        [cat_x, cat_y, loc_x, loc_y, logq_x_c, logq_y_c] ...
            = SchooseCatsCoupled(state_x, state_y, i, newage_x, newage_y);
    else
        [cat_x, cat_y, loc_x, loc_y] = deal([]);
        [logq_x_c, logq_y_c] = deal(0);
    end

    logq_x = logq_x_t + logq_x_c;
    logq_y = logq_y_t + logq_y_c;
end

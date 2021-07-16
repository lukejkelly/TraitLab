function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = MoveCatLocCoupled(state_x, state_y)
    % Sample catastrophes to move and new locations from maximal coupling
    if state_x.ncat > 0 && state_y.ncat > 0
        [i_x, i_y, j_x, j_y] = sampleCatIndexCoupled(state_x, state_y);
        catloc_new = rand;
        [nstate_x, U_x, OK_x, logq_x] ...
            = MoveCatLoc.getOutputsNew(state_x, i_x, j_x, catloc_new);
        [nstate_y, U_y, OK_y, logq_y] ...
            = MoveCatLoc.getOutputsNew(state_y, i_y, j_y, catloc_new);
    else
        % Let marginal move deal with it
        [nstate_x, U_x, OK_x, logq_x] = MoveCatLoc(state_x);
        [nstate_y, U_y, OK_y, logq_y] = MoveCatLoc(state_y);
    end
end

function [nstate, U, OK, logq] = MoveCatLoc(state)
    % Move catastrophe location j along branch i
    if state.ncat > 0
        [i, j] = sampleCatIndex(state);
        catloc_new = rand;
        [nstate, U, OK, logq] ...
            = MoveCatLoc.getOutputsNew(state, i, j, catloc_new);
    else
        [nstate, U, OK, logq] = MoveCatLoc.getOutputsFail(state);
    end
end

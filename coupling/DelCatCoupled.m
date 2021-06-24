function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = DelCatCoupled(state_x, state_y)
    % Maximally coupled version of core/AddCat; see its comments for details
    global BORROWING
    if state_x.ncat == 0 || state_y.ncat == 0
        % Let core/DelCat deal with it
        [nstate_x, U_x, OK_x, logq_x] = DelCat(state_x);
        [nstate_y, U_y, OK_y, logq_y] = DelCat(state_y);
    else
        % Attempt maximally coupled deletion
        if BORROWING
            % Catastrophes individually identifiable
            [i_x, i_y, j_x, j_y] = sampleCatIndexCoupled(state_x, state_y);
        else
            [i_x, i_y] = sampleBranchProportionalToCatCountCoupled(state_x, ...
                                                                   state_y);
            [j_x, j_y] = deal([]);
        end
        [nstate_x, U_x, OK_x, logq_x] = getOutputs(state_x, i_x, j_x);
        [nstate_y, U_y, OK_y, logq_y] = getOutputs(state_y, i_y, j_y);
    end
end

function [nstate, U, OK, logq] = getOutputs(state, i, j)
    global BORROWING
    nstate = state;
    nstate.cat(i) = state.cat(i) - 1;
    nstate.ncat = state.ncat - 1;

    if BORROWING
        nstate.tree(i).catloc(j) = [];
    end

    U = above(i, state.tree, state.root);
    OK = 1;

    if BORROWING
        logq = log(state.ncat) - log(state.length);
    else
        dt = nstate.tree(nstate.tree(i).parent).time - nstate.tree(i).time;
        logq = -log(state.length) + log(dt) + log(state.ncat) - log(state.cat(i));
    end
end

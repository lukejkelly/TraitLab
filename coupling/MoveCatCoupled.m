function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = MoveCatCoupled(state_x, state_y)
    % Maximally coupled version of core/MoveCat
    global BORROWING

    if state_x.ncat == 0 || state_y.ncat == 0
        % Let core/MoveCat deal with it
        [nstate_x, U_x, OK_x, logq_x] = MoveCat(state_x);
        [nstate_y, U_y, OK_y, logq_y] = MoveCat(state_y);
    else
        % Attempt maximally coupled deletion
        [old_x, old_y] = sampleBranchByCatCountCoupled(state_x, state_y);
        if BORROWING
            [j_x, j_y] = sampleCatIndexCoupled(state_x, state_y, old_x, old_y);
            catloc_new = rand;
        else
            [j_x, j_y, catloc_new] = deal([]);
        end

        root = state_x.root;
        [new_x, new_y, q1_x, q1_y, q2_x, q2_y] = GetLegalCoupled(...
            state_x.tree, state_y.tree, old_x, old_y, root);

        [nstate_x, U_x, OK_x, logq_x] ...
            = getOutputs(state_x, old_x, new_x, q1_x, q2_x, j_x, catloc_new);
        [nstate_y, U_y, OK_y, logq_y] ...
            = getOutputs(state_y, old_y, new_y, q1_y, q2_y, j_y, catloc_new);
    end
end

function [nstate, U, OK, logq] = getOutputs(state, old, new, q1, q2, j, catloc_new)
    global BORROWING
    nstate = state;
    nstate.cat(old) = state.cat(old) - 1;
    nstate.cat(new) = state.cat(new) + 1;
    if BORROWING
        nstate.tree(old).catloc(j) = [];
        nstate.tree(new).catloc = sort([nstate.tree(new).catloc, catloc_new]);
    end

    U = above([old, new], nstate.tree, nstate.root);
    OK = 1;

    if BORROWING
        dt1 = nstate.tree(nstate.tree(old).parent).time - nstate.tree(old).time;
        dt2 = nstate.tree(nstate.tree(new).parent).time - nstate.tree(new).time;
        logq = log(q1) - log(q2) - log(dt1) + log(dt2);
    else
        logq = log(q1) - log(q2) - log(state.cat(old)) + log(nstate.cat(new));
    end
end

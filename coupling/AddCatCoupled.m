function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = AddCatCoupled(state_x, state_y)
    % Maximally coupled version of core/AddCat; see comments there for details
    global BORROWING

    [r_x, r_y] = sampleBranchProportionalToLengthCoupled(state_x, state_y);

    if BORROWING
        new_catloc = rand;
    else
        new_catloc = [];
    end

    [nstate_x, U_x, OK_x, logq_x] = getOutputs(state_x, r_x, new_catloc);
    [nstate_y, U_y, OK_y, logq_y] = getOutputs(state_y, r_y, new_catloc);
end

function [nstate, U, OK, logq] = getOutputs(state, r, new_catloc)
    global BORROWING
    nstate = state;
    nstate.cat(r) = state.cat(r) + 1;
    nstate.ncat = state.ncat + 1;
    U = above(r, nstate.tree, nstate.root);
    OK = 1;
    if BORROWING
        nstate.tree(r).catloc = sort([nstate.tree(r).catloc, new_catloc]);
        logq = log(state.length) - log(nstate.ncat);
    else
        dt = state.tree(state.tree(r).parent).time - state.tree(r).time;
        logq = log(state.length) - log(dt) - log(nstate.ncat) ...
               + log(nstate.cat(r));
    end
end

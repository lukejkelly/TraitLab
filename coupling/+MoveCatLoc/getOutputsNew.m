function [nstate, U, OK, logq] = getOutputsNew(state, i, j, catloc_new)
    % Replace catastrophe location and get associated quantities
    catloc = state.tree(i).catloc;
    catloc(j) = catloc_new;

    nstate = state;
    nstate.tree(i).catloc = sort(catloc);

    U = above(i, nstate.tree, nstate.root);
    OK = 1;
    logq = 0;
end

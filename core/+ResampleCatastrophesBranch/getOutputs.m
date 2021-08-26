function [nstate, U, OK, logq] = getOutputs(state, i, cat, catloc)
    % Create output variables after resampling catastrophes
    global BORROWING

    nstate = state;
    nstate.cat(i) = cat;
    nstate.ncat = sum(nstate.cat);

    if BORROWING && cat > 0
        % Facilitates marginal and coupled approaches
        nstate.tree(i).catloc = sort(catloc(1:cat));
    else
        nstate.tree(i).catloc = [];
    end

    U = above(i, nstate.tree, nstate.root);
    OK = 1;

    logq = ResampleCatastrophesBranch.logPrior(state, i) ...
           - ResampleCatastrophesBranch.logPrior(nstate, i);
end

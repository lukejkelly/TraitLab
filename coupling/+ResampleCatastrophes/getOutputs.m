function [nstate, U, OK, logq] = getOutputs(state, cat, catloc)
    % Create output variables after resampling catastrophes
    global BORROWING

    nstate = state;
    nstate.cat(:) = cat;
    nstate.ncat = sum(cat);

    if BORROWING
        for j = 1:(2 * nstate.NS)
            if nstate.cat(j) > 0
                % Facilitates marginal and coupled approaches
                nstate.tree(j).catloc = sort(catloc{j}(1:nstate.cat(j)));
            else
                nstate.tree(j).catloc = [];
            end
        end
    end

    U = above(nstate.leaves, nstate.tree, nstate.root);
    OK = 1;

    logq = ResampleCatastrophes.logPrior(state) ...
           - ResampleCatastrophes.logPrior(nstate);
end

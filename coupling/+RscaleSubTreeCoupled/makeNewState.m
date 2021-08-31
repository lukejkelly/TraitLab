function [nstate, U, logq, OK] = makeNewState(state, i, lf, nd, t0, newage)
    % Based on core/RscaleSubTree
    nstate = state;
    OK = 1;

    % Update internal node times
    rho = (newage - t0) / (state.tree(i).time - t0);
    for j = nd
       nstate.tree(j).time = t0 + rho * (nstate.tree(j).time - t0);
    end

    % Fail if any leaf or i older than its parent
    for j = lf
       if nstate.tree(nstate.tree(j).parent).time < nstate.tree(j).time
          OK = 0;
          break;
       end
    end
    k = nstate.tree(i).parent;
    if nstate.tree(i).time > nstate.tree(k).time
        OK = 0;
    end

    if OK
       U = above(nd, nstate.tree, nstate.root);
       logq = (length(nd) - 2) * log(rho);
       nstate.length = TreeLength(nstate.tree, nstate.root);
    else
       U = [];
       logq = 0;
    end
end

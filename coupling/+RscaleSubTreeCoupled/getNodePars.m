function [lf, nd, t0] = getNodePars(state, nu)
    % Get subtree leaf and internal node sets and youngest leaf
    global LEAF;
    Vel = nu(1, :);
    lf = [];
    nd = [];
    for j = Vel
       if state.tree(j).type == LEAF
          lf = [lf, j];
       else
          nd = [nd, j];
       end
    end
    t0 = min([state.tree(lf).time]);
end

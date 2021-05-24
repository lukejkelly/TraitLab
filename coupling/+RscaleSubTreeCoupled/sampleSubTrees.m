function [i_x, i_y, nu_x, nu_y] = sampleSubTrees(state_x, state_y)
    % Sample subtree roots, also return list of nodes beneath them and parents
    global LEAF;

    r = state_x.root;
    if r ~= state_y.root
        error('Root indices do not match');
    end

    nu_x = progeny(state_x.tree, r, LEAF);
    nu_y = progeny(state_y.tree, r, LEAF);

    rp = @() RscaleSubTreeCoupled.rDiscrete(nu_x);
    ldp = @(i) RscaleSubTreeCoupled.lpDiscrete(i, nu_x);

    rq = @() RscaleSubTreeCoupled.rDiscrete(nu_y);
    ldq = @(i) RscaleSubTreeCoupled.lpDiscrete(i, nu_y);

    [i_x, i_y] = maximalCouplingLog(rp, ldp, rq, ldq);

    if any([state_x.tree(i_x).type, state_y.tree(i_y).type] == LEAF)
        error('Cannot have leaf as subtree root');
    end

    if i_x ~= r
       nu_x = progeny(state_x.tree, i_x, LEAF);
    end
    if i_y ~= r
       nu_y = progeny(state_y.tree, i_y, LEAF);
    end
end

function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledXi(state_x, state_y, mcmc)

    % Previous version for testing against
    nstate_x = state_x;
    nstate_y = state_y;

    % Same leaf indices in both x and y after housekeeping
    leaf = randsample(state_x.leaves, 1);
    % state_x.leaves(ceil(rand * length(state.leaves)));

    % xi'_tilde = 1 - xi' = U(del, del + deldel) * (1 - xi)
    [nstate_x_tree_leaf_xi_tilde, nstate_y_tree_leaf_xi_tilde] ...
        = maximalCouplingUniformScaling(...
            1 - state_x.tree(leaf).xi, ...
            1 - state_y.tree(leaf).xi, ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    nstate_x.tree(leaf).xi = 1 - nstate_x_tree_leaf_xi_tilde;
    nstate_y.tree(leaf).xi = 1 - nstate_y_tree_leaf_xi_tilde;

    var_x = (1 - nstate_x.tree(leaf).xi) / (1 - state_x.tree(leaf).xi);
    var_y = (1 - nstate_y.tree(leaf).xi) / (1 - state_y.tree(leaf).xi);

    logq_x = -log(var_x);
    logq_y = -log(var_y);

    U_x = above(leaf, state_x.tree, state_x.root);
    U_y = above(leaf, state_y.tree, state_y.root);

    OK_x = nstate_x.tree(leaf).xi >= 0;
    OK_y = nstate_y.tree(leaf).xi >= 0;

end

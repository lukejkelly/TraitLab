function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledXi(state_x, state_y, mcmc)

    % Same leaf indices in both x and y after housekeeping
    leaf = randsample(state_x.leaves, 1);
    % state_x.leaves(ceil(rand * length(state.leaves)));

    % xi'_tilde = 1 - xi' = U(del, del + deldel) * (1 - xi)
    [xi_tilde_x, xi_tilde_y] ...
        = maximalCouplingUniformScaling(...
            1 - state_x.tree(leaf).xi, ...
            1 - state_y.tree(leaf).xi, ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    [nstate_x, U_x, OK_x, logq_x] = getOutputs(state_x, leaf, xi_tilde_x);
    [nstate_y, U_y, OK_y, logq_y] = getOutputs(state_y, leaf, xi_tilde_y);

end

function [nstate, U, OK, logq] = getOutputs(state, leaf, xi_tilde)
    nstate = state;
    nstate.tree(leaf).xi = 1 - xi_tilde;
    U = above(leaf, state.tree, state.root);
    OK = nstate.tree(leaf).xi >= 0;
    var = (1 - nstate.tree(leaf).xi) / (1 - state.tree(leaf).xi);
    logq = -log(var);
end

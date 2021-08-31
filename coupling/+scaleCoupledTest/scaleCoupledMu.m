function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledMu(state_x, state_y, mcmc)

    % Previous version for testing against
    global DEPNU DONTMOVECATS;
    nstate_x = state_x;
    nstate_y = state_y;

    [nstate_x.mu, nstate_y.mu] ...
        = maximalCouplingUniformScaling(...
            state_x.mu, ...
            state_y.mu, ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    var_x = nstate_x.mu / state_x.mu;
    var_y = nstate_y.mu / state_y.mu;

    logq_x = -log(var_x);
    logq_y = -log(var_y);

    U_x = state_x.nodes;
    U_y = state_y.nodes;

    if DEPNU
        nstate_x.nu = state_x.nu / var_x;
        nstate_y.nu = state_y.nu / var_y;
    end

    OK_x = ~(DONTMOVECATS && nstate_x.mu < 1e-5);
    OK_y = ~(DONTMOVECATS && nstate_y.mu < 1e-5);

end

function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledKappa(state_x, state_y, mcmc)

    % Previous version for testing against
    global DEPNU;
    nstate_x = state_x;
    nstate_y = state_y;

    [nstate_x.kappa, nstate_y.kappa] ...
        = maximalCouplingUniformScaling(...
            state_x.kappa, ...
            state_y.kappa, ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    var_x = nstate_x.kappa / state_x.kappa;
    var_y = nstate_y.kappa / state_y.kappa;

    logq_x = -log(var_x);
    logq_y = -log(var_y);

    U_x = state_x.nodes;
    U_y = state_y.nodes;

    OK_x = ((0.25 <= nstate_x.kappa) && (nstate_x.kappa <= 1));
    OK_y = ((0.25 <= nstate_y.kappa) && (nstate_y.kappa <= 1));

    if OK_x && DEPNU
        nstate_x.nu = state_x.nu * var_x;
    end
    if OK_y && DEPNU
        nstate_y.nu = state_y.nu * var_y;
    end

end

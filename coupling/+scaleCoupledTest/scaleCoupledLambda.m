function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledLambda(state_x, state_y, mcmc)

    % Previous version for testing against
    global DEPNU;
    nstate_x = state_x;
    nstate_y = state_y;
    [OK_x, OK_y] = deal(1);

    [nstate_x.lambda, nstate_y.lambda] ...
        = maximalCouplingUniformScaling(...
            state_x.lambda, ...
            state_y.lambda, ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    var_x = nstate_x.lambda / state_x.lambda;
    var_y = nstate_y.lambda / state_y.lambda;

    logq_x = -log(var_x);
    logq_y = -log(var_y);

    U_x = state_x.nodes;
    U_y = state_y.nodes;

    if DEPNU
        nstate_x.nu = state_x.nu * var_x;
        nstate_y.nu = state_y.nu * var_y;
    end

end

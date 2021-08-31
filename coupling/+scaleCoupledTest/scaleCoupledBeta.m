function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledBeta(state_x, state_y, mcmc)

    % Previous version for testing against
    nstate_x = state_x;
    nstate_y = state_y;
    [OK_x, OK_y] = deal(1);

    [nstate_x.beta, nstate_y.beta] ...
        = maximalCouplingUniformScaling(...
            state_x.beta, ...
            state_y.beta, ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    var_x = nstate_x.beta / state_x.beta;
    var_y = nstate_y.beta / state_y.beta;

    logq_x = -log(var_x);
    logq_y = -log(var_y);

    U_x = state_x.nodes;
    U_y = state_y.nodes;

end

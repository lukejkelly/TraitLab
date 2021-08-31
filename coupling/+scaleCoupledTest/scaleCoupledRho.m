function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledRho(state_x, state_y, mcmc)

    % Previous version for testing against
    nstate_x = state_x;
    nstate_y = state_y;
    [OK_x, OK_y] = deal(1);

    [nstate_x.rho, nstate_y.rho] ...
        = maximalCouplingUniformScaling(...
            state_x.rho, ...
            state_y.rho, ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    var_x = nstate_x.rho / state_x.rho;
    var_y = nstate_y.rho / state_y.rho;

    logq_x = -log(var_x);
    logq_y = -log(var_y);

    U_x = [];
    U_y = [];

end

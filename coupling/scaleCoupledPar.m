function [nstate_x, nstate_y, logq_x, logq_y, var_x, var_y] ...
        = scaleCoupledPar(state_x, state_y, mcmc, par)
    % Sample from coupling of state_*.par and return associated quantities
    nstate_x = state_x;
    nstate_y = state_y;

    [nstate_x.(par), nstate_y.(par)] ...
        = maximalCouplingUniformScaling(...
            state_x.(par), ...
            state_y.(par), ...
            mcmc.update.del, ...
            mcmc.update.del + mcmc.update.deldel);

    [logq_x, var_x] = getOutputs(state_x, nstate_x, par);
    [logq_y, var_y] = getOutputs(state_y, nstate_y, par);
end

function [logq, var] = getOutputs(state, nstate, par)
    var = nstate.(par) / state.(par);
    logq = -log(var);
end

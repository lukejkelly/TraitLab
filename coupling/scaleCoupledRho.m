function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledRho(state_x, state_y, mcmc)

    [nstate_x, nstate_y, logq_x, logq_y, ~, ~] ...
        = scaleCoupledPar(state_x, state_y, mcmc, 'rho');

    [U_x, OK_x] = getOutputs();
    [U_y, OK_y] = getOutputs();

end

function [U, OK] = getOutputs()
    U = [];
    OK = 1;
end

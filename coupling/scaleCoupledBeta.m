function [nstate_x, nstate_y, U_x, U_y, OK_x, OK_y, logq_x, logq_y] ...
        = scaleCoupledBeta(state_x, state_y, mcmc)

    [nstate_x, nstate_y, logq_x, logq_y, ~, ~] ...
        = scaleCoupledPar(state_x, state_y, mcmc, 'beta');

    [U_x, OK_x] = getOutputs(nstate_x);
    [U_y, OK_y] = getOutputs(nstate_y);

end

function [U, OK] = getOutputs(nstate)
    U = nstate.nodes;
    OK = 1;
end

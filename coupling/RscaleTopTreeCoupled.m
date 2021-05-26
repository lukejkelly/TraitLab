function [nstate_x, nstate_y, U_x, U_y, logq_x, logq_y, OK_x, OK_y] ...
        = RscaleTopTreeCoupled(state_x, state_y, prior, del, deldel)

    % RscaleTopTree always includes root so reuse RscaleCoupled to get var terms
    [var_x, var_y] = RscaleCoupled(state_x, state_y, del, del + deldel);

    % Apply updates separately to X and Y replacing del by var_ and deldel by 0
    % to erase randomness in RscaleTopTree
    [nstate_x, U_x, ~, logq_x, OK_x] = RscaleTopTree(state_x, prior, var_x, 0);
    [nstate_y, U_y, ~, logq_y, OK_y] = RscaleTopTree(state_y, prior, var_y, 0);

end

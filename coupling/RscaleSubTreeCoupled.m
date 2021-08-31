function [nstate_x, nstate_y, U_x, U_y, logq_x, logq_y, OK_x, OK_y] ...
        = RscaleSubTreeCoupled(state_x, state_y, del, deldel)
    % Sample subtree roots and get associated parameters
    [i_x, i_y, nu_x, nu_y] = RscaleSubTreeCoupled.sampleSubTrees(...
        state_x, state_y);

    % Subtree leaf and node sets and youngest leaf
    [lf_x, nd_x, t0_x] = RscaleSubTreeCoupled.getNodePars(state_x, nu_x);
    [lf_y, nd_y, t0_y] = RscaleSubTreeCoupled.getNodePars(state_y, nu_y);

    % t' ~ U(a, b) = U(t0 + del * (t - t0), t0 + (del + deldel) * (t - t0))
    [a_p, b_p] = RscaleSubTreeCoupled.getAgePars(state_x, i_x, t0_x, del, deldel);
    [a_q, b_q] = RscaleSubTreeCoupled.getAgePars(state_y, i_y, t0_y, del, deldel);

    [newage_x, newage_y] = maximalCouplingUniform(a_p, b_p, a_q, b_q);

    % Update subtree times
    [nstate_x, U_x, logq_x, OK_x] = RscaleSubTreeCoupled.makeNewState(...
        state_x, i_x, lf_x, nd_x, t0_x, newage_x);
    [nstate_y, U_y, logq_y, OK_y] = RscaleSubTreeCoupled.makeNewState(...
        state_y, i_y, lf_y, nd_y, t0_y, newage_y);
end

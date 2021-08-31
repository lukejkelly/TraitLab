function [ncat_x, ncat_y, cat_x, cat_y, loc_x, loc_y, logq_x, logq_y] = BcatsCoupled(...
        state_x, state_y, i, j_x, j_y, k_x, k_y, newage_x, newage_y)
    % How we modify catastrophes in an SPR move depends on whether pa(i), j  or
    % neither is currently the root,
    %     — x will be the new catastrophe count for either pa(i) or sib(i) in
    %       the proposed state
    %     - logq is the log-ratio of proposals to account for this, the terms
    %       corresponding to locations on branches are handled by
    %       catastropheScalingFactor
    % i's parent p is to be regrafted onto edge from j to parent k at newage
    % Currently, i's sibling is h and p's parent is q
    global MCMCCAT
    if ~MCMCCAT || any(isempty([newage_x, newage_y]))
        error('MCMCCAT = %i, any(isempty([newage_x, newage_y])) = %i', ...
              MCMCCAT, any(isempty([newage_x, newage_y])));
    end

    [rf_x, lf_x, lb_x] = Bcats.getDistributions(state_x, i, j_x, k_x, newage_x);
    [rf_y, lf_y, lb_y] = Bcats.getDistributions(state_y, i, j_y, k_y, newage_y);

    [jn_x, jn_y] = maximalCouplingLog(rf_x, lf_x, rf_y, lf_y);

    logq_x = lb_x - lf_x(jn_x);
    logq_y = lb_y - lf_y(jn_y);

    [ncat_x, cat_x] = Bcats.newCatCounts(state_x, i, j_x, jn_x);
    [ncat_y, cat_y] = Bcats.newCatCounts(state_y, i, j_y, jn_y);

    cat = BcatsCoupled.getOverallCounts(cat_x, cat_y);
    loc = Bcats.newCatLocations(cat);
    loc_x = BcatsCoupled.getSpecificLocations(cat_x, loc);
    loc_y = BcatsCoupled.getSpecificLocations(cat_y, loc);
end

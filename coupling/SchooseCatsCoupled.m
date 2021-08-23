function [cat_x, cat_y, loc_x, loc_y, logq_x, logq_y] = SchooseCatsCoupled(...
        state_x, state_y, i, newage_x, newage_y)
    % Branch i the same in both and due to housekeeping root index is the same
    global MCMCCAT
    if ~MCMCCAT || any(isempty([newage_x, newage_y]))
        error('MCMCCAT = %i, any(isempty([newage_x, newage_y])) = %i', ...
              MCMCCAT, any(isempty([newage_x, newage_y])));
    end

    [catc_x, ncat_x, wc_x, wn_x] = SchooseCats.branchesAndWeights(...
        state_x, i, newage_x);
    [catc_y, ncat_y, wc_y, wn_y] = SchooseCats.branchesAndWeights(...
        state_y, i, newage_y);

    [rf_x, lf_x, lb_x] = SchooseCats.getDistributions(...
        catc_x, ncat_x, wc_x, wn_x);
    [rf_y, lf_y, lb_y] = SchooseCats.getDistributions(...
        catc_y, ncat_y, wc_y, wn_y);

    % We only attempt to couple if catastrophe counts are identical as mnpdf
    % gets number of trials from summing argument
    if ncat_x == ncat_y
        [catn_x, catn_y] = maximalCouplingLog(rf_x, lf_x, rf_y, lf_y);
    else
        catn_x = rf_x();
        catn_y = rf_y();
    end
    cat_x = SchooseCats.catOutput(catn_x);
    cat_y = SchooseCats.catOutput(catn_y);

    % Similar to BcatsCoupled: get overall catastrophe counts and couple
    % sampling locations (if necessary)
    cat = SchooseCatsCoupled.getOverallCounts(cat_x, cat_y);
    loc = Bcats.newCatLocations(cat);
    loc_x = SchooseCatsCoupled.getSpecificLocations(cat_x, loc);
    loc_y = SchooseCatsCoupled.getSpecificLocations(cat_y, loc);

    % Contribution from locations dealt with by catastropheScalingFactor
    logq_x = lb_x - lf_x(catn_x);
    logq_y = lb_y - lf_y(catn_y);
end

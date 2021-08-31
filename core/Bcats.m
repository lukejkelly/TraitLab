function [ncat, cat, loc, logq] = Bcats(state, i, j, k, newage)
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
    if ~MCMCCAT || isempty(newage)
        error('MCMCCAT = %i, isempty(newage) = %i', MCMCCAT, isempty(newage));
    end
    [rf, lf, lb] = Bcats.getDistributions(state, i, j, k, newage);
    jn = rf();
    logq = lb - lf(jn);

    [ncat, cat] = Bcats.newCatCounts(state, i, j, jn);
    loc = Bcats.newCatLocations(cat);
end

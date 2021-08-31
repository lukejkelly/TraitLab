function [cat, loc, logq] = SchooseCats(state, i, newage)
    % How we modify catastrophes in when resampling a node time: Multinomial
    % shuffling of catastrophes with weights proportional to branch lengths
    global MCMCCAT
    if ~MCMCCAT || isempty(newage)
        error('MCMCCAT = %i, isempty(newage) = %i', MCMCCAT, isempty(newage));
    end

    [catc, ncat, wc, wn] = SchooseCats.branchesAndWeights(state, i, newage);

    [rf, lf, lb] = SchooseCats.getDistributions(catc, ncat, wc, wn);
    catn = rf();

    cat = SchooseCats.catOutput(catn);
    loc = Bcats.newCatLocations(cat);
    logq = lb - lf(catn);
end

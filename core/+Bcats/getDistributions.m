function [rf, lf, lb] = getDistributions(state, i, j, k, newage)
    % handles to simulate and evaluate log-probability for forward move, and
    % log-probability for return move (whatever they may be)
    % return move does not depend on proposed state
    Root = state.root;

    % handles to sample new cat count on j (and possibly p)
    [jt, jc] = getTimeAndCats(state, j);
    [kt, ~] = getTimeAndCats(state, k);

    dj_new = newage - jt;
    if j == Root
        % p will become new root and lose catastrophes so forward move is
        % sample from prior on branch from j to new parent p
        rf = @() Bcats.samplePriorCounts(dj_new, state.rho);
        lf = @(z) Bcats.logPriorCounts(z, dj_new, state.rho);
    else
        % forward move is binomial split of j's cat count to get  (jc', pc')
        wf = dj_new / (kt - jt);
        rf = @() binornd(jc, wf);
        lf = @(z) log(binopdf(z, jc, wf));
    end

    % log-probability of proposal to return to original state
    [p, q, h] = Bcats.pqh(state, i);

    [oldage, pc] = getTimeAndCats(state, p);
    [qt, ~] = getTimeAndCats(state, q);
    [ht, hc] = getTimeAndCats(state, h);

    dh_old = oldage - ht;
    if p == Root
        % backwards move is sample from prior to restore h's cat count
        lb = Bcats.logPriorCounts(hc, dh_old, state.rho);
    else
        % backwards move is binomial to restore h's and p's cat counts
        wb = dh_old / (qt - ht);
        lb = log(binopdf(hc, hc + pc, wb));
    end
end

function [it, ic] = getTimeAndCats(state, i)
    it = state.tree(i).time;
    ic = state.cat(i);
end

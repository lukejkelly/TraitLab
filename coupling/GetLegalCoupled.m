function [new_x, new_y, q1_x, q1_y, q2_x, q2_y] = GetLegalCoupled(s_x, s_y, ...
        old_x, old_y, root)
    % Maximally  coupled version of core/GetLegal
    [poss_x, q1_x] = GetLegalCoupled.getPoss(s_x, old_x, root);
    [poss_y, q1_y] = GetLegalCoupled.getPoss(s_y, old_y, root);

    [rp, ldp] = getDistributionTerms(poss_x, q1_x);
    [rq, ldq] = getDistributionTerms(poss_y, q1_y);
    [new_x, new_y] = maximalCouplingLog(rp, ldp, rq, ldq);

    [~, q2_x] = GetLegalCoupled.getPoss(s_x, new_x, root);
    [~, q2_y] = GetLegalCoupled.getPoss(s_y, new_y, root);
end

function [r, ld] = getDistributionTerms(poss, q)
    r = @() GetLegalCoupled.rNeighbour(poss, q);
    ld = @(x) GetLegalCoupled.ldNeighbour(x, poss, q);
end

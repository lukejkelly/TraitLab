function [new_minage, kT, logq] = getUniformParameters(i, j, k, s, THETA)
    global OTHER ROOT

    iT = s(i).time;
    iP = s(i).parent;

    jT = s(j).time;
    kT = s(k).time;

    new_minage = max(iT, jT);
    new_range = kT - new_minage;

    CiP = s(iP).child(OTHER(s(i).sibling));
    CiPT = s(CiP).time;

    % Sampling is the same but logq depends on whether iP was root
    if s(iP).type == ROOT
        logq = (CiPT - s(iP).time) * THETA + log(new_range) + log(THETA);
    else
        PiP = s(iP).parent;
        PiPT = s(PiP).time;
        old_minage = max(iT, CiPT);
        old_range = PiPT - old_minage;

        logq = log(new_range) - log(old_range);
    end
end

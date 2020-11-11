function [logq] = getCoupling1Parameter(i, jT, s, newage, THETA)
    global OTHER

    delta = newage - jT;

    iT = s(i).time;
    iP = s(i).parent;

    PiP = s(iP).parent;
    CiP = s(iP).child(OTHER(s(i).sibling));

    PiPT = s(PiP).time;
    CiPT = s(CiP).time;

    old_minage = max(iT, CiPT);
    old_range = PiPT - old_minage;

    logq = delta * THETA - log(THETA * old_range);
end

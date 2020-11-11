function [newage, logq] = sampleMarginal(i, j, k, s, THETA)
    global OTHER ROOT

    iT = s(i).time;
    iP = s(i).parent;

    if s(j).type == ROOT

       jT = s(j).time;
       delta = -(1 / THETA) * log(rand);
       newage = jT + delta;

       PiP = s(iP).parent;
       CiP = s(iP).child(OTHER(s(i).sibling));

       PiPT = s(PiP).time;
       CiPT = s(CiP).time;
       old_minage = max(iT,CiPT);
       old_range = PiPT - old_minage;

       q = exp(delta * THETA) / (THETA * old_range);
       logq = log(q);

    elseif s(iP).type == ROOT

       jT = s(j).time;
       kT = s(k).time;
       new_minage = max(iT, jT);
       new_range = kT - new_minage;
       newage = new_minage + rand * new_range;

       CiP = s(iP).child(OTHER(s(i).sibling));
       CiPT = s(CiP).time;
       q = exp((CiPT - s(iP).time) * THETA) * new_range * THETA;
       logq = log(q);

    else

       jT = s(j).time;
       kT = s(k).time;
       new_minage = max(iT, jT);
       new_range = kT - new_minage;
       newage = new_minage + rand * new_range;

       PiP = s(iP).parent;
       CiP = s(iP).child(OTHER(s(i).sibling));

       PiPT = s(PiP).time;
       CiPT = s(CiP).time;
       old_minage = max(iT, CiPT);
       old_range = PiPT - old_minage;

       q = new_range / old_range;
       logq = log(q);
    end
end

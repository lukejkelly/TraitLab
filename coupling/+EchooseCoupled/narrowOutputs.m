function [j, iP, jP, OK] = narrowOutputs(i, state, prior)
    % Taken from corresponding part of core/Echoose
    global OTHER

    s = state.tree;
    iP = s(i).parent;
    jP = s(iP).parent;
    j = s(jP).child(OTHER(s(iP).sibling));

    if (s(j).time < s(iP).time) ...
            && (~prior.isclade || (prior.isclade && isequal(s(iP).clade, ...
                                                            s(jP).clade)))
        OK = 1;
    else
        OK = 0;
    end
end

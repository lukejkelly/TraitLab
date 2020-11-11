function [j, k, FAIL] = getNarrowDestination(i, s)
    % <k, j> = <pa(pa(i)), sib(pa(i))>
    global ROOT OTHER
    iP = s(i).parent;
    if s(iP).type == ROOT
        FAIL = 1;
        k = -1;
        j = -1;
    else
        FAIL = 0;
        k = s(iP).parent;
        j = s(k).child(OTHER(s(iP).sibling));
    end
end

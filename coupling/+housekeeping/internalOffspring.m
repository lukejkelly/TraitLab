function r = internalOffspring(s, i)
    % Internal nodes strictly beneath i in tree s
    global LEAF ANST
    p = progeny(s, i, LEAF);
    q = setdiff(p(1, :), i);
    r = q([s(q).type] == ANST);
end

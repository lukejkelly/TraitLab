function nstate = updateCats(state, i, j, ncat, cat, loc)
    % Update catastrophe counts and locations in advance of SPR topology move
    global BORROWING
    nstate = state;
    nstate.ncat = ncat;

    [p, ~, h] = Bcats.pqh(state, i);

    nstate.cat(j) = cat.j;
    nstate.cat(h) = cat.h;
    nstate.cat(p) = cat.p;

    if BORROWING
        nstate.tree(j).catloc = sort(loc.j);
        nstate.tree(h).catloc = sort(loc.h);
        nstate.tree(p).catloc = sort(loc.p);
    end
end

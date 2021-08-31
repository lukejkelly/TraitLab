function state = moveCats(state, i, cat, loc)
    % Update catastrophe counts (and locations) after move on node time
    global BORROWING

    b = [i, state.tree(i).child];
    state.cat(b(1)) = cat.i;
    state.cat(b(2)) = cat.j;
    state.cat(b(3)) = cat.k;

    if BORROWING
        state.tree(b(1)).catloc = sort(loc.i);
        state.tree(b(2)).catloc = sort(loc.j);
        state.tree(b(3)).catloc = sort(loc.k);
    end
end

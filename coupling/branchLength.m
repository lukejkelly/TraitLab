function bl = branchLength(state, i)
    % Branch length if below root, 0 otherwise
    global ROOT
    if state.tree(i).type < ROOT
        bl = state.tree(state.tree(i).parent).time - state.tree(i).time;
    else
        bl = 0;
    end
end

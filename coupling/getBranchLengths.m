function bl = getBranchLengths(state)
    % Returns a vector where each entry is the length of the branch from that
    % node to its parent but ignoring the root and Adam branches
    global ROOT
    L = 2 * state.NS;
    bl = zeros(1, L);
    for i = 1:L
        if state.tree(i).type < ROOT
            bl(i) = state.tree(state.tree(i).parent).time - state.tree(i).time;
        end
    end
end

function j = getIndex(state, i, k)
    j = find(ismembertol(state.tree(i).catloc, k));
    if length(j) > 1
        error('Only one catastrophe location should be selected');
    end
end

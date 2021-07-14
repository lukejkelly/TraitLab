function [i, j] = getIndex(state, k)
    i = find(cellfun(@(x) ismembertol(k, x), {state.tree.catloc}));
    if length(i) > 1
        error('Only one catastrophe branch should be selected');
    end
    j = find(ismembertol(state.tree(i).catloc, k));
    if length(j) > 1
        error('Only one catastrophe location should be selected');
    end
end

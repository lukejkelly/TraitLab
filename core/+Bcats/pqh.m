function [p, q, h] = pqh(state, i)
    % Return indices of i's parent, grandparent and sibling
    global OTHER
    s = state.tree;
    p = s(i).parent;
    q = s(p).parent;
    h = s(p).child(OTHER(s(i).sibling));
end

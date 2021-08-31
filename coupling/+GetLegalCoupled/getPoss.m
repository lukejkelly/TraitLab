function [poss, q] = getPoss(s, i, root)
    global OTHER
    poss = setdiff([s(i).child, s(i).parent, ...
                    s(s(i).parent).child(OTHER(s(i).sibling))], ...
                   root);
    q = length(poss);
    if q == 0
        error('Node does not have a sibling nor children and its parent is the root.');
    end
end

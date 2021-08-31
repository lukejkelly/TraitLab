function ld = ldNeighbour(x, poss, q)
    if ismember(x, poss)
        ld = -log(q);
    else
        ld = -Inf;
    end
end

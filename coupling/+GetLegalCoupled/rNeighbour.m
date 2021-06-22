function new = rNeighbour(poss, q)
    if q ~= length(poss)
        % TODO: This error should never arise, remove once debugging over
        error('Length of poss and q do not match');
    end
    new = poss(ceil(rand * q));
end

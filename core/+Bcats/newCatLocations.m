function loc = newCatLocations(cat)
    % Sample catastrophe locations on branches given counts
    global BORROWING
    if BORROWING
        loc = structfun(@(n) rand(n > 0, n), cat, 'UniformOutput', false);
    else
        loc = [];
    end
end

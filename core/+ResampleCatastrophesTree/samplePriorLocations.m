function catloc = samplePriorLocations(cat)
    % Sample catastrophe locations on branches given counts
    global BORROWING
    if BORROWING
        catloc = arrayfun(@(n) rand(n > 0, n), cat, 'UniformOutput', false);
    else
        catloc = [];
    end
end

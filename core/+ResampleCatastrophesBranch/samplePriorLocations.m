function catloc = samplePriorLocations(cat)
    % Sample catastrophe locations on branches given counts
    global BORROWING
    if BORROWING && cat > 0
        catloc = rand(1, cat);
    else
        catloc = [];
    end
end

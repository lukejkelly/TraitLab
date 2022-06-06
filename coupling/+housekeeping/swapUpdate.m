function [t, j] = swapUpdate(t, j, i1, i2)
    % Swap nodes i1 and i2 in tree t and the corresponding entries in vector j
    t = housekeeping.swapNodes(t, i1, i2);
    j = housekeeping.swapEntries(j, i1, i2);
end

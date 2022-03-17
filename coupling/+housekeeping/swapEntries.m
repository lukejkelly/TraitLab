function j = swapEntries(j, i1, i2)
    % Swap entries i1 and i2 in vector j
    j([i1, i2]) = j([i2, i1]);
end

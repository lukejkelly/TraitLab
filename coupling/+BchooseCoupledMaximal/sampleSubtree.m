function i = sampleSubtree(s)
    global ROOT
    N = length(s) - 1;
    i = ceil(N * rand);
    while s(i).type == ROOT
       i = ceil(N * rand);
    end
end

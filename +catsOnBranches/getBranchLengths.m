function d = getBranchLengths(outPath, outFile)
    global ROOT
    tAll = catsOnBranches.getTrees(outPath, outFile);
    t = rnextree(tAll{1});
    d = zeros(1, length(t));
    for j = 1:size(d, 2)
        if t(j).type < ROOT
            d(j) = t(t(j).parent).time - t(j).time;
        end
    end
end

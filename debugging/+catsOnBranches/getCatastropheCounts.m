function n = getCatastropheCounts(outPath, outFile, burnin)
    global ROOT
    cAll = catsOnBranches.getTrees(outPath, [outFile, 'cat']);
    n = zeros(length(cAll) - burnin, length(rnextree(cAll{1})));
    for i = 1:size(n, 1)
        c = rnextree(cAll{i + burnin});
        for j = 1:size(n, 2)
            if c(j).type < ROOT
                n(i, j) = round(c(c(j).parent).time - c(j).time - 0.1);
            end
        end
    end
end

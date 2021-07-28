function d = getBranchLengths(outPath, outFile, burnin)
    tAll = catsOnBranches.getTrees(outPath, outFile);
    t = rnextree(tAll{1});
    n = length(tAll) - burnin;
    m = length(t);
    d = zeros(n, m);
    for i = 1:n
        state = struct('NS', m / 2, tree = rnextree(tAll{i + burnin}));
        d(i, :) = getBranchLengths(state);
    end
end

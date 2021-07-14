function tAll = getTrees(outPath, outFile)
    tOut = readrunfile(sprintf('%s/%s.nex', outPath, outFile));
    tInds = cellfun(@(x) ~isempty(x), tOut{2});
    tAll = tOut{2}(tInds);
end

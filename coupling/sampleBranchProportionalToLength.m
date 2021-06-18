function r = sampleBranchProportionalToLength(state)
    bl = getBranchLengths(state);
    r = find(rand < cumsum(bl) ./ sum(bl), 1);
end

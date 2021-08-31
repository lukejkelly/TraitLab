function r = sampleBranchByCatCount(state)
    if sum(state.cat) ~= state.ncat
        error('Inconsistent catastrophe counts');
    end
    r = find(rand < cumsum(state.cat) ./ state.ncat, 1);
end

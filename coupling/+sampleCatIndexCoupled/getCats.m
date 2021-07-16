function [locs, inds] = getCats(state)
    % Used for unit testing index-sampling functions
    % locs is a vector of catastrophe (relative) locations, the first column of
    % inds is the corresponding catastrophe branches and the second column the
    % location indices
    locs = nan(state.ncat, 1);
    inds = nan(state.ncat, 2);
    k = 1;
    for j = find(state.cat(:)')
        jInds = k:(k + state.cat(j) - 1);
        locs(jInds) = state.tree(j).catloc;
        inds(jInds, 1) = j;
        inds(jInds, 2) = 1:state.cat(j);
        k = k + state.cat(j);
    end
end

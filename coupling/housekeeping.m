% We assume that s1 was only ever a valid modification of s2 so leaves and Adam
% have the same indices in both, as do all nodes within a clade
% TODO: We should implement more general housekeeping without these assumptions

function nstate2 = housekeeping(state1, state2)
    global LEAF ANST ROOT
    % Match indices of nodes with common descendent leaves (indexed by .Name)
    s1 = state1.tree;
    s2 = state2.tree;
    t2 = s2;
    r = {s1([s1.type] == 0).Name};  % Reference list of leaf names
    a1 = housekeeping.getLeafArray(s1, r);
    i1 = (1:length(s1))';
    done = false;
    while ~done
        a2 = housekeeping.getLeafArray(t2, r);
        [c1, c2] = ismember(a1, a2, 'rows');
        i2 = find(c1 & (i1 ~= c2), 1);  % i2 = c1 & (i1 ~= c2);
        if any(i2)
            % i = find(i2, 1);
            % t2:c2(i2) has same role as s1:i2 so we swap t2:i2 and t2:c2(i2)
            t2 = housekeeping.swapNodes(t2, i2, c2(i2));
        else
            done = true;
        end
    end
    % Pairs of common siblings have same ordering
    t2 = housekeeping.matchSiblings(s1, t2);
    % TODO: remove this after further testing
    if equaltrees(s2, t2) == 0
        save('coupling/housekeeping-state.mat');
        error('Trees do not match after housekeeping');
    end
    % Update state variables
    nstate2 = state2;
    nstate2.tree = t2;
    if all(sort(find([t2.type] == LEAF)) == sort(state1.leaves))
        % TODO: check if ordering matters
        nstate2.leaves = state1.leaves;
    else
        error('Leaf sets do not match after housekeeping');
    end
    if find([t2.type] == ROOT) == state1.root
        % TODO: check if ordering matters
        nstate2.root = state1.root;
    else
        error('Roots do not match after housekeeping');
    end
    if all(sort([find([t2.type] == ANST), nstate2.root]) == sort(state1.nodes))
        % TODO: check if ordering matters
        nstate2.nodes = state1.nodes;
    else
        error('Internal nodes do not match after housekeeping');
    end
    nstate2.cat(:) = cellfun(@length, {nstate2.tree.catloc});
    nstate2 = UpdateClades(nstate2, [nstate2.leaves, nstate2.nodes], ...
                           length(nstate2.claderoot));
    if ~all(isequaln(nstate2.claderoot, state1.claderoot))
        save('coupling/housekeeping-error.mat')
        error('Clade roots do not match after housekeeping');
    end
end

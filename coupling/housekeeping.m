function nstate2 = housekeeping(state1, state2)
    % Match indices of nodes with common descendent leaves (indexed by .Name)
    % For each tree, we identify which leaves (ordered by reference list r) are
    % beneath each node using getLeafArray then permute node indices with
    % swapNodes
    % Common subtrees have the same indices in both trees
    % getLeafArray does not consider the Adam node as if its index is out of
    % position then it will get swapped anyway
    % We then update sibling information with matchSiblings, if possible, as
    % well as any clades or catastrophes

    global LEAF ANST ROOT BORROWING

    s1 = state1.tree;
    s2 = state2.tree;
    t2 = s2;

    % Permute indices so that subtrees have common parent indices
    r = {s1([s1.type] == 0).Name};
    a1 = housekeeping.getLeafArray(s1, r);
    i1 = (1:length(s1))';
    j2 = 1:length(s2);
    done = false;
    while ~done
        a2 = housekeeping.getLeafArray(t2, r);
        [c1, c2] = ismember(a1, a2, 'rows');
        i2 = find(c1 & (i1 ~= c2), 1);
        if any(i2)
            % This move also swaps non-relationship information such as Name and
            % xi; that is, t(j).Name is now s(k).Name. We do not use internal
            % node Name or xi information anywhere (that I know of) and leaves
            % should never swap with internal nodes in MCMC (which uses a
            % similar operation and has been doing so for years) or housekeeping
            % so I think we can safely ignore it for checking coupling
            t2 = housekeeping.swapNodes(t2, i2, c2(i2));
            j2([i2, c2(i2)]) = j2([c2(i2), i2]);
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

    nstate2.cat = nstate2.cat(j2);
    if BORROWING && any(nstate2.cat ~= cellfun(@length, {nstate2.tree.catloc}'))
        error('Catastrophes do not match');
    end

    nstate2 = UpdateClades(nstate2, [nstate2.leaves, nstate2.nodes], ...
                           length(nstate2.claderoot));
    if ~all(isequaln(nstate2.claderoot, state1.claderoot))
        save('coupling/housekeeping-error.mat')
        error('Clade roots do not match after housekeeping');
    end

     % Update node likelihood information
    if ~BORROWING
        nstate2 = MarkRcurs(nstate2, [nstate2.leaves, nstate2.nodes], 1);
    end
    % if BORROWING && ~ismembertol(nstate2.loglkd, logLkd2(nstate2))
    %     keyboard;
    % elseif ~BORROWING && ~ismembertol(nstate2.loglkd, LogLkd(nstate2))
    %     keyboard;
    % elseif check(state1) || check(state2) || check(nstate2)
    %     keyboard;
    % end
end

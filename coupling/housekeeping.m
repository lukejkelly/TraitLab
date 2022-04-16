function nstate2 = housekeeping(state1, state2)
    % Match indices of nodes with common descendent leaves (indexed by .Name)
    % We identify which leaves (ordered by reference list r) are beneath each
    % node using getLeafArray then permute node indices with swapNodes, finally
    % update sibling information with matchSiblings, if necessary, as well as
    % any clades or catastrophes
    %  -- Common clades use the same indices in both trees
    %  -- Indices of leaves and Adam should be unchanged throughout MCMC
    %  -- getLeafArray does not consider the Adam node as if its index is out of
    %     position then it will get swapped anyway

    global LEAF ANST ROOT BORROWING

    s1 = state1.tree;
    s2 = state2.tree;
    t2 = s2;

    % Permute indices so that subtrees have common parent indices
    r = {s1([s1.type] == LEAF).Name};
    a1 = housekeeping.getLeafArray(s1, r, state1.root);
    i1 = (1:length(s1))';
    j2 = 1:length(s2);
    done = false;
    while ~done
        root2 = find([t2.type] == ROOT);
        a2 = housekeeping.getLeafArray(t2, r, root2);
        [c1, c2] = ismember(a1, a2, 'rows');
        i2 = find(c1 & (i1 ~= c2), 1);
        if any(i2)
            [t2, j2] = housekeeping.swapUpdate(t2, j2, i2, c2(i2));
        else
            done = true;
        end
    end

    % Common ancestral nodes
    for k = state1.nodes(c1(state1.nodes))
        done = false;
        while ~done
            k1 = housekeeping.internalOffspring(s1, k);
            k2 = housekeeping.internalOffspring(t2, k);
            [l, l1, l2] = setxor(k1, k2);
            if any(l)
                [t2, j2] = housekeeping.swapUpdate(t2, j2, k1(l1(1)), k2(l2(1)));
            else
                done = true;
            end
        end
    end

    % Pairs of common siblings have same ordering
    t2 = housekeeping.matchSiblings(s1, t2);
    % % TODO: remove this after further testing
    % if equaltrees(s2, t2) == 0
    %     save('coupling/housekeeping-state.mat');
    %     error('Trees do not match after housekeeping');
    % end

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

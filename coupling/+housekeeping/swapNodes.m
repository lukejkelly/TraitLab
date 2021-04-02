function t = swapNodes(s, j, k)
    % Swap subtree-parent role of nodes j and k
    % ns_j and ns_k are the new sibling roles for j and k after swapping
    % Almost identical to core/swap.m but preserving sibling order

    % TODO: remove after further testing
    if j == k
        error('We should not be trying to swap a node with itself');
    end

    % Order so only special case is j has k as parent
    if isequaln(s(k).parent, j)
        % If k's parent is j then swap j and k
        [j, k] = deal(k, j);
    end

    sj = s(j);
    sk = s(k);

    % Exchange node fields
    t = s;
    t(j) = sk;
    t(k) = sj;

    % Update nodes depending on how j and k relate
    if isempty(sj.parent) || sj.parent ~= k
        % k is not j's parent in s
        if ~isempty(sj.parent)
            t(sj.parent).child(sj.sibling) = k;
        end
        if ~isempty(sj.child)
            [t(sj.child).parent] = deal(k);
        end
        if ~isempty(sk.parent)
            t(sk.parent).child(sk.sibling) = j;
        end
        if ~isempty(sk.child)
            [t(sk.child).parent] = deal(j);
        end
    else
        % k is j's parent in s
        if ~isempty(sj.child)
            [t(sj.child).parent] = deal(k);
        end
        t(k).parent = j;
        t(j).parent = sk.parent;
        if ~isempty(sk.parent)
            t(sk.parent).child(sk.sibling) = j;
            t(j).sibling = sk.sibling;
        end
        t(j).child(sj.sibling) = k;
        c = sk.child(sk.child ~= j);
        if ~isempty(c)
            t(c).parent = j;
        end
    end
    % TODO: remove this after further testing
    result = checktree(t, size(s, 2) / 2);
    if ~isempty(result)
        save('coupling/housekeeping-error.mat')
        error('Not a valid tree');
    end
end

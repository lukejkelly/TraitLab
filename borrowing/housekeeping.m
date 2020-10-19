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
    a1 = getLeafArray(s1, r);
    i1 = (1:length(s1))';
    done = false;
    while ~done
        a2 = getLeafArray(t2, r);
        [c1, c2] = ismember(a1, a2, 'rows');
        i2 = find(c1 & (i1 ~= c2), 1);  % i2 = c1 & (i1 ~= c2);
        if any(i2)
            % i = find(i2, 1);
            % t2:c2(i2) has same role as s1:i2 so we swap t2:i2 and t2:c2(i2)
            t2 = swapNodes(t2, i2, c2(i2));
        else
            done = true;
        end
    end
    % TODO: remove this after further testing
    if equaltrees(s2, t2) == 0
        disp('Trees do not match after housekeeping');
        keyboard;
    end
    % Update state variables
    nstate2 = state2;
    nstate2.tree = t2;
    nstate2.leaves = find([t2.type] == LEAF);
    nstate2.root = find([t2.type] == ROOT);
    nstate2.nodes = [find([t2.type] == ANST), nstate2.root];
    nstate2.cat(:) = cellfun(@length, {nstate2.tree.catloc});
    nstate2 = UpdateClades(nstate2, [nstate2.leaves, nstate2.nodes], ...
                           length(nstate2.claderoot));
    if ~all(isequaln(nstate2.claderoot, state1.claderoot))
        disp('Clade roots do not match after housekeeping');
        keyboard;
    end
end

function a = getLeafArray(s, r, a, l)
    % Starting at node l and working towards the leaves of s, a indicates the
    % the leaves (identified by name) beneath each node according to the
    % reference list r
    global LEAF ROOT;

    if nargin == 2
        a = zeros(length(s), length(r));
        l = find([s.type] == ROOT);
    end

    if s(l).type == LEAF
        a(l, ismember(r, s(l).Name)) = 1;
    else
        a1 = getLeafArray(s, r, a, s(l).child(1));
        a2 = getLeafArray(s, r, a, s(l).child(2));
        a(~all(a1 == 0, 2), :) = a1(~all(a1 == 0, 2), :);
        a(~all(a2 == 0, 2), :) = a2(~all(a2 == 0, 2), :);
        a(l, :) = sum(a(s(l).child, :));
    end
end

function t = swapNodes(s, j, k)
    % Swap subtree-parent role of nodes j and k
    % ns_j and ns_k are the new sibling roles for j and k after swapping
    % Almost identical to core/swap.m but preserving sibling order

    % Only special case is j has k as parent
    if isequaln(s(k).parent, j)
        [j, k] = deal(k, j);
    end

    t = s;
    sj = s(j);
    sk = s(k);

    % Exchange node fields
    t(j) = sk;
    t(k) = sj;

    % Update nodes depending on how j and k relate
    if ~isequaln(sj.parent, k)
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
        disp('Not a valid tree');
        keyboard;
    end
end

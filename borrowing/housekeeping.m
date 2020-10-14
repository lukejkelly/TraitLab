function [t2] = housekeeping(s1, s2)
    % Match indices of nodes with common descendent leaves (indexed by .Name)
    t2 = s2;
    r = {s1([s1.type] == 0).Name}; % {s1.Name};
    a1 = getLeafArray(s1, r);
    i1 = (1:length(s1))';
    done = false;
    while ~done
        a2 = getLeafArray(t2, r);
        [c1, c2] = ismember(a1, a2, 'rows');
        i2 = c1 & (i1 ~= c2);
        if any(i2)
            i = find(i2, 1);
            % t2:c2(i) is playing same role as s1:i so t2:i and t2:c2(i)
            t2 = swapNodes(t2, i, c2(i));
        else
            done = true;
        end
    end
end

function [a] = getLeafArray(s, r, a, l)
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
    result = checktree(t, size(s, 2) / 2);
    if ~isempty(result)
        disp('error pp');
        keyboard;
        pause;
    end
end

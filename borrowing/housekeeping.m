function [t2] = housekeeping(s1, s2)
    % Match indices of nodes with common descendent leaves (indexed by .Name)
    t2 = s2;
    r = {s1.Name};
    a1 = getLeafArray(s1, r);
    i1 = (1:length(s1))';
    done = 0;
    while ~done
        a2 = getLeafArray(t2, r);
        [c1, c2] = ismember(a1, a2, 'rows');
        i2 = c1 & (i1 ~= c2);
        if any(i2)
            i = find(i2, 1);
            % We have c2(i) playing i's role, but want i playing i's role
            fprintf('switch %i and %i\n', i, c2(i));
            t2 = swapNodes(t2, i, c2(i));
        else
            done = 1;
        end
    end
end

% Exchange nodes j and k in tree s
function t = swapNodes(s, j, k)
    % Swap subtree-parent role of nodes i and j in when doing housekeeping
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
        c = sk.child(sk.child ~= j);
        % t(j).child = [k, c];
        % t(k).sibling = 1;
        t(j).child(sj.sibling) = k;
        if ~isempty(c)
            % t(c).sibling = 2;
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

function [a] = getLeafArray(s, r, a, l)
    % Starting at node l and working towards the leaves of s, a indicates the
    % the leaves beneath the node according to the reference list r
    % each node
    % Leaves identified by their name
    % TODO: identify leaves by their data
    global LEAF ROOT;

    if nargin == 2
        a = zeros(length(s));
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

% function [a] = getLeafArray(s, a, l)
%     % Starting at node l and working towards the leaves, list the leaves beneath
%     % each node
%     % Leaves identified by their index
%     global LEAF ROOT;
%
%     if nargin == 1
%         a = zeros(length(s));
%         l = find([s.type] == ROOT);
%     end
%
%     if s(l).type == LEAF
%         a(l, l) = 1;
%     else
%         a1 = getLeafArray(s, a, s(l).child(1));
%         a2 = getLeafArray(s, a, s(l).child(2));
%         a(~all(a1 == 0, 2), :) = a1(~all(a1 == 0, 2), :);
%         a(~all(a2 == 0, 2), :) = a2(~all(a2 == 0, 2), :);
%         a(l, :) = sum(a(s(l).child, :));
%     end
% end

function a = getLeafArray(s, r, a, l)
    % Starting at node l (root or lower) and working towards the leaves of s,
    % a indicates the leaves (identified by name) beneath each node according to
    % the reference list r; that is, a is a |s|x|a| array with a_ij = 1 if leaf
    % r_j is beneath node s_i

    global LEAF ROOT;

    if nargin == 2
        a = zeros(length(s), length(r));
        l = find([s.type] == ROOT);
    end

    if s(l).type == LEAF
        a(l, ismember(r, s(l).Name)) = 1;
    else
        a1 = housekeeping.getLeafArray(s, r, a, s(l).child(1));
        a2 = housekeeping.getLeafArray(s, r, a, s(l).child(2));
        a(~all(a1 == 0, 2), :) = a1(~all(a1 == 0, 2), :);
        a(~all(a2 == 0, 2), :) = a2(~all(a2 == 0, 2), :);
        a(l, :) = sum(a(s(l).child, :));
    end
end

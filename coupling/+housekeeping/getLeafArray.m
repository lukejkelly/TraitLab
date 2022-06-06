function a = getLeafArray(s, r, l)
    % Starting at node l (root or lower) and working towards the leaves of s,
    % a indicates the leaves (identified by name) beneath each node according to
    % the reference list r; that is, a is a |s|x|a| array with a_ij = 1 if leaf
    % r_j is beneath node s_i

    if isempty(s(l).child)
        a = zeros(length(s), length(r));
        a(l, ismember(r, s(l).Name)) = 1;
    else
        a1 = housekeeping.getLeafArray(s, r, s(l).child(1));
        a2 = housekeeping.getLeafArray(s, r, s(l).child(2));
        a = a1 + a2;
        a(l, :) = sum(a(s(l).child, :));
    end
end

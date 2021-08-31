function [iT, kT, jT, a, b] = nodeTimesAndRanges(i, state)
    % Times of node i, its parent and its oldest child, ranges for sampling
    % new root time
    s = state.tree;
    iT = s(i).time;
    kT = s(s(i).parent).time;
    jT = max(s(s(i).child).time);
    a = (iT + jT) / 2;
    b = 2 * iT - jT;
end

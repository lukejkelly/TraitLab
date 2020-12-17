function s = sampleSyntheticTree(L, root_time)
    % Node times scaled so that root age is root_time
    global ROOT ANST

    THETA = 1;
    s = ExpTree(L, THETA);

    old_root_time = s([s.type] == ROOT).time;
    for j = find([s.type] == ANST | [s.type] == ROOT)
       s(j).time = s(j).time * root_time / old_root_time;
    end
end

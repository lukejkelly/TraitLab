function v = valid(i, r, s)
    iT = s(i).time;
    j = r(:)';
    k = [s(j).parent];
    kT = [s(k).time];
    v = j(~(kT <= iT | i == j | i == k));
end

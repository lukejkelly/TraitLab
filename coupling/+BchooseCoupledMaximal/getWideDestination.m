function [j, k, FAIL] = getWideDestination(i, r, N, s)
    if N > 4
        iT = s(i).time;
        rInd = 1;
        j = r(rInd);
        k = s(j).parent;
        while (s(k).time <= iT || i == j || i == k)
            rInd = rInd + 1;
            j = r(rInd);
            k = s(j).parent;
        end
        FAIL = any([j, k] == s(i).parent);
    else
        j = -1;
        k = -1;
        FAIL = 1;
    end
end

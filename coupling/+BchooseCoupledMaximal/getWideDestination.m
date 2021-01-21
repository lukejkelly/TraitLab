function [j, k, FAIL] = getWideDestination(i, r, N, s)
    % From a vector of possible destination edges (randomised in th calling
    % function) return the first valid
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
        if any([j, k] == s(i).parent)
            FAIL = 1;
        else
            FAIL = 0;
        end
    else
        j = -1;
        k = -1;
        FAIL = 1;
    end
end

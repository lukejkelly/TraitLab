function [x, y] = maximalCouplingLog(rp, ldp, rq, ldq)
    % Sampling from a maximal coupling of p and q
    % Sampling functions rp and rq, ldp and ldq are corresponding log densities
    x = rp();
    if log(rand) + ldp(x) <= ldq(x)
        y = x;
    else
        y = nan;
        while isnan(y)
            ys = rq();
            if log(rand) + ldq(ys) > ldp(ys)
                y = ys;
            end
        end
    end
end

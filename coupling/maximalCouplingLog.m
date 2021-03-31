function [x, y] = maximalCouplingLog(rp, ldp, rq, ldq)
    % Sampling from a maximal coupling of x ~ p and y ~ q
    % rp and rq generate samples, ldp and ldq are log densities
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

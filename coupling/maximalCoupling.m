function [x, y] = maximalCoupling(rp, dp, rq, dq)
    % Sampling from a maximal coupling of p and q
    % rp generates a draw from p, dp is the corresponding density; q likewise
    x = rp();
    if rand * dp(x) <= dq(x)
        y = x;
    else
        y = nan;
        n = 0;
        while isnan(y)
            ys = rq();
            if rand * dq(ys) > dp(ys)
                y = ys;
            end
            n = n + 1;
            if n == 1000
                keyboard;
            end
        end
    end
end

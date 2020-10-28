function [x, y] = maximalCouplingLog(rp, ldp, rq, ldq)
    % Sampling from a maximal coupling of p and q
    % Densities ldp and ldq are on log scale
    % TODO: test
    x = rp();
    if log(rand) + ldp(x) <= ldq(x)
        y = x;
    else
        y = nan;
        % n = 0;
        while isnan(y)
            ys = rq();
            if log(rand) + ldq(ys) > ldp(ys)
                y = ys;
            end
            % n = n + 1;
            % if n == 1000
            %     keyboard;
            % end
        end
    end
end

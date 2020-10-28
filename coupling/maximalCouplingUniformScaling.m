function [x, y] = maximalCouplingUniformScaling(xc, yc, a, b)
    % Sample from a maximal coupling of p = U(a*xc, b*xc) and q = U(a*xc, b*xc)

    % Sampling and log density functions
    rp = @() prototypeSample(xc, a, b);
    rq = @() prototypeSample(yc, a, b);

    ldp = @(x) prototypeLogDensity(x, xc, a, b);
    ldq = @(y) prototypeLogDensity(y, yc, a, b);

    % Avoid while loop and sampling from piecewise
    if ismembertol(xc, yc)
        % Already coupled
        [x, y] = deal(rp());
    elseif ~(a * xc < b * yc || a * yc < b * xc)
        % No overlap
        x = rp();
        y = rq();
    elseif xc < yc
        % Overlap and p > q on supp p
        x = rp();
        % If x on supp min(p, q) then attempt to couple, otherwise sample y
        % uniformly on supp q \ supp p
        if log(rand) + ldp(x) <= ldq(x)
            y = x;
        else
            % x not on supp p and q so sample y uniformly on supp q \ supp p
            y = b * (xc + rand * (yc - xc));
        end
    else
        % Overlap and q > p on supp q
        y = rq();
        % If y on supp min(p, q) then attempt to couple, otherwise sample x
        % uniformly on supp p \ supp q
        if log(rand) + ldq(y) <= ldp(y)
            x = y;
        else
            x = b * (yc + rand * (xc - yc));
        end
    end
end


% Prototypes for coupling
function u = prototypeSample(v, a, b)
    % Sample u ~ Unif(va, vb)
    u = v * (a + rand * (b - a));
end

% function d = prototypeDensity(u, v, a, b)
%     % Density of u ~ Unif(va, vb)
%     if (v * a <= u && u <= v * b)
%         d = 1 / (v * (b - a));
%     else
%         d = 0;
%     end
% end

function ld = prototypeLogDensity(u, v, a, b)
    % Log density of u ~ Unif(va, vb)
    if (v * a <= u && u <= v * b)
        ld = -(log(v) + log(b - a));
    else
        ld = -Inf;
    end
end

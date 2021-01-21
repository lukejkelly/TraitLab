function [x, y] = maximalCouplingUniformShiftedExponential(a, b, c, theta)
    % Sampling from a maximal coupling of p = U(a, b) and q = ShiftExp(c, theta)
    rp = @() sampleUniform(a, b);
    rq = @() sampleShiftedExponential(c, theta);

    ldp = @(x) logDensityUniform(x, a, b);
    ldq = @(y) logDensityShiftedExponential(y, c, theta);

    [x, y] = maximalCouplingLog(rp, ldp, rq, ldq);
end

function u = sampleUniform(a, b)
    u = a + rand * (b - a);
end

function ld = logDensityUniform(u, a, b)
    if (a <= u && u <= b)
        ld = -log(b - a);
    else
        ld = -Inf;
    end
end

function u = sampleShiftedExponential(c, theta)
    u = c - log(rand) / theta;
end

function ld = logDensityShiftedExponential(u, c, theta)
    if u < c
        ld = -Inf;
    else
        ld = log(theta) - theta * (u - c);
    end
end

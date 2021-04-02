function [x, y] = maximalCouplingUniform(a_p, b_p, a_q, b_q)
    % Sampling from a maximal coupling of p = U(a_p, b_p) and q = U(a_q, b_q)

    % Sampling and log density functions
    rp = @() prototypeSample(a_p, b_p);
    rq = @() prototypeSample(a_q, b_q);

    ldp = @(x) prototypeLogDensity(x, a_p, b_p);
    ldq = @(y) prototypeLogDensity(y, a_q, b_q);

    [x, y] = maximalCouplingLog(rp, ldp, rq, ldq);
end

% Prototypes for coupling
function u = prototypeSample(a, b)
    % Sample u ~ Unif(a, b)
    u = a + rand * (b - a);
end

function ld = prototypeLogDensity(u, a, b)
    % Log density of u ~ Unif(a, b)
    if (a <= u && u <= b)
        ld = -log(b - a);
    else
        ld = -Inf;
    end
end

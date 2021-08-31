function [x, y] = maximalCouplingShiftedExponential(a_p, a_q, theta)
    % Sample (x, y) from a maximal coupling of p(x - a_p) = Exp(theta) and ...
    % q(y -  a_q) = Exp(theta)

    % Overlap
    log_min_PQ = -theta * abs(a_p - a_q);

    if log(rand) <= log_min_PQ
        % Coupling
        x = sampleShiftedExponential(max(a_p, a_q), theta);
        y = x;
    else
        % Sample independently
        if a_p < a_q
            x = sampleShiftedTruncatedExponential(a_p, a_q, theta);
            y = sampleShiftedExponential(a_q, theta);
        else
            x = sampleShiftedExponential(a_p, theta);
            y = sampleShiftedTruncatedExponential(a_q, a_p, theta);
        end
    end
end

function t = sampleShiftedExponential(a, l)
    % t = a + Exp(l)
    t = a - log(rand) / l;
end

function t = sampleShiftedTruncatedExponential(a, b, l)
    % t = a + Exp(l) on (a, b)
    t = a - log(1 - rand * (1 - exp(-l * (b - a)))) / l;
end

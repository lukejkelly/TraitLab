function i = rDiscrete(nu)
    % Adapted from core/RscaleSubTree
    vec = cumsum(nu(2, :));
    tot = vec(end);
    r = rand * tot;
    i = nu(1, find(r < vec, 1));
end

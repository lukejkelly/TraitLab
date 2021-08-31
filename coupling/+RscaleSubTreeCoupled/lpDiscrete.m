function lp = lpDiscrete(i, nu)
    % Adapted from core/RscaleSubTree
    num = nu(2, nu(1, :) == i);
    den = sum(nu(2, :));
    lp = log(num) - log(den);
end

function lp = lpdfG(x, a, b)
% Evaluate log PDF at point x of Gamma with shape a and rate b

if any([x, a, b] <= 0)
    error('arguments must be non-negative: x = %g, a = %g, b = %g', x, a, b);
end

lp = a * log(b) + (a - 1) * log(x) - b * x - gammaln(a);

end

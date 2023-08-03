function f = randG(a, b)
% Generates a Gamma(shape = a, rate = b) random variable f with PDF
%     p(f) = b^a f^(a - 1) * exp(-b * f) * / Gamma(a)
% by Marsaglia and Tsang (2000) https://dl.acm.org/doi/10.1145/358407.358414
% RJR, 12/02/2009

if (a <= 0)
    error('first parameter must be non-negative: shape a = %g', a);
elseif a < 1
    f = randG(a + 1, b) * rand^(1 / a);
elseif (a == 1)
    f = -log(rand) / b;
else % a > 1
	d = a - 1 / 3;
	c = 1 / sqrt(9 * d);

	while 1
		x = randn;
		v = (1 + c * x)^3;
		if (v > 0) && (log(rand) < 0.5 * x^2 + d - d * v + d * log(v))
            break;
        end
	end
	f = d * v / b;
end

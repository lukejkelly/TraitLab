function lp = logpost(lambda, state, x)
% evaluate Gamma(shape = n, rate = x) marginal log posterior at lambda
%     n: the number of recorded traits
%     x: the sum of the expected pattern frequencies for registered traits
% prior on lambda is improper 1 / lambda and lambda scales expected pattern
% frequencies in Poisson likelihood

if state.L == 0
    warning('n = 0: improper prior on lambda');
    lp = -log(lambda);
else
    lp = lpdfG(lambda, state.L, x);
end

end

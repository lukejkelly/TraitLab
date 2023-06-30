function lambda = sample(n, x)
% sample from Gamma(n, x) marginal posterior for lambda, where
%     n: the number of recorded traits
%     x: the sum of the expected pattern frequencies for registered traits
% prior on lambda is improper 1 / lambda and lambda scales expected pattern
% frequencies in Poisson likelihood

if n == 0
    warning('cannot sample from improper prior on lambda, return 0 instead');
    lambda = 0;
else
    lambda = randG(n, x);
end

end

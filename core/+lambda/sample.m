function lambda = sample(state, x)
% sample from Gamma(n, x) marginal posterior for lambda, where
%     n: the number of recorded traits
%     x: the sum of the expected pattern frequencies for registered traits
% prior on lambda is improper 1 / lambda and lambda scales expected pattern
% frequencies in Poisson likelihood

if state.L == 0
    lambda = 0;
else
    lambda = randG(state.L, x);
end

end

function [lp, k, theta] = LogRhoPrior(rho)
% returns log-density lp on prior rho ~ Gamma(shape = k, scale = theta)
% note that the parameterisation to sample a Gamma RV via randG or evaluate its
% log PDF via lpdfG is shape = k and rate = 1 / theta
% prior 5% and 95% quantiles are [0.000035, 0.00078] so corresponding quantiles
% for time between catastrophes is [1300, 28000] years

k = 1.5;
theta = 0.0002;

lp = lpdfG(rho, k, 1 / theta);

end

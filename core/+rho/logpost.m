function lp = logpost(rho, c, d)
% evaluate Gamma(shape = a + c, rate = b + d) marginal posterior at rho, where
%     c: number of catastrophes on the tree
%     d: length of tree below root
% a and b are parameters of Gamma(shape = a, rate = b) prior on rho

[~, a, theta] = LogRhoPrior(0);
b = 1 / theta;
lp = lpdfG(rho, a + c, b + d);

end

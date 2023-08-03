function lp = logpost(rho, state)
% evaluate Gamma(shape = a + c, rate = b + d) marginal posterior at rho, where
%     a, b: parameters of Gamma(shape = a, rate = b) prior
%     c: number of catastrophes on the tree
%     d: length of tree below root

[~, a, theta] = LogRhoPrior(0);
b = 1 / theta;
c = state.ncat;
d = state.length;

lp = lpdfG(rho, a + c, b + d);

end

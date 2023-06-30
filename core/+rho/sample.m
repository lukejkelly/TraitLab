function rho = sample(c, d)
% sample from Gamma(shape = a + c, rate = b + d) marginal posterior on rho
%     c: number of catastrophes on the tree
%     d: length of tree below root
% a and b are parameters of Gamma(shape = a, rate = b) prior

[~, a, theta] = LogRhoPrior(0);
b = 1 / theta;
rho = randG(a + c, b + d);

end

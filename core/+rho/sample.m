function rho = sample(state)
% draw rho from its Gamma(shape = a + c, rate = b + d) marginal posterior
%     a, b: parameters of Gamma(shape = a, rate = b) prior
%     c: number of catastrophes on the tree
%     d: length of tree below root

[~, a, theta] = LogRhoPrior(0);
b = 1 / theta;
c = state.ncat;
d = state.length;

rho = randG(a + c, b + d);

end

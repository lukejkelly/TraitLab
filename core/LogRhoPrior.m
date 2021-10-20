function [p k theta]=LogRhoPrior(rho)

% p=LogRhoPrior(rho) returns the log-prior of rho
% notice that the parameterisation to sample from randG is
% randG(k,1/theta), also the 5,95 % quantiles for k=1.5, theta=0.0002 
% are at about [0.000035,0.00078] ao that's one catastrophe every [1300,28000].
% k=1.5; theta=0.0002; N=10000; v=1:N; for i=v, v(i)=randG(k,1/theta); end
% u=sort(v); 1./u(round([0.05,0.95].*N)),1/mean(u)

%GKN Feb 08 - changed RJR's original
k=1.5;
theta=.0002;

p=(k-1)*log(rho)-rho/theta-k*log(theta)-gammaln(k);
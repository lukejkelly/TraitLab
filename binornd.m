function R=binornd(n,p,k)

%%function R=binornd(n,p,k)
%GKN implementation of a Binomial sampler

if nargin==2, k=1; end
if (nargin>3 | length(k)>1), error('this isnt the MatLab Stats binornd()'); end

R=sum(rand(n,k)<p);






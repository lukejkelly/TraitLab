function x = randln(mu,sigma,ismeanandvar,m,n)
% function x = randln(mu,sigma,ismeanandvar,m,n)
% generate m by n matrix of lognormal random variables
% if ismeanandvar = 1, then X ~ LN(a,b) such that E(X)=mu and
% var(X)=sigma^2
% if ismeanandvar = 0, X ~ LN(mu,sigma)
% m,n are optional inputs, assumed 1 if not given

if nargin == 3
    m = 1; n = 1;
elseif nargin == 4
    n = 1;
elseif nargin ~= 5
    error('Number of arguments in randln must be 3, 4 or 5')
end

if sigma < 0
    error('sigma must be >= 0 in randln')
end

if ismeanandvar
    if mu <= 0
        error('mu must be greater than 0 in randln when ismeanandvar = 1')
    end
    b = log(1+(sigma/mu)^2);
    a = log(mu) - b/2;
    
else
    a = mu;
    b = sigma;
end

x = exp(sqrt(b)*randn(m,n) + a);
function p=LogMuPrior(mu,ALPHA,BETA)

%p=zeros(size(mu));
%p(mu<ALPHA | mu>BETA)=-inf;

p=-log(mu);
p(mu<ALPHA | mu>BETA)=-inf;

%p=(ALPHA-1).*log(mu)-mu./BETA;

%unnormalised gamma
%mean is alpha*beta var is alpha*beta^2

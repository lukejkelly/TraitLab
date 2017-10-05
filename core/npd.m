function f=npd(prior,ActI,mu,p,t)

global YULE FLAT

a=exp(-mu.*t);

X=(2*(1-a)+1-(1-p.*a).^2)./mu;
Y=1-a+p.*a.*(1-p.*a);
Z=(p.*a).^2;

L=length([ActI{:}]);

nDif=length(ActI{1})+length(ActI{2});
nSam=length(ActI{3});

f = -L*log(mu) - log(L) - L.*log(X) + nDif.*log(Y) + nSam.*log(Z);

if prior.type==YULE
   f=f-log(2.*t);
elseif ~(prior.type==FLAT)
   disp('PRIOR type not recognised in pd.m');
   keyboard;pause;
end

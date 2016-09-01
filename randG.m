function f=randG(a,b)
% Generates a Gamma random variable with pdf
% x^(a-1)*exp(-c*x)*c*a/Gamma(a)
% Following method by Marsaglia and Tsang (2000)
% RJR, 12/02/2009

if (a<=0)

   disp('error in randG at point 1'); 

elseif a<1
	f=randG(a+1,b)*rand^(1/a);

elseif (a==1)
  f=-log(rand)/b;

else % a>1
	d=a-1/3;
	c=1/(sqrt(9*d));

	while 1
		x=randn;
		v=(1+c*x)^3;
		if (v>0) && (log(rand)<.5*x^2+d-d*v+d*log(v))
            break;
        end
	end;

	f=d*v/b;
end

function UpdateHEt(mu)

BACmu=mu;
load IEexample;
mu=BACmu;

for i=1:(NS-1)
   for j=(i+1):NS
      HEt(i,j)=NaiveMLE(s(i).ActI{3},s(j).ActI{3},mu);
   end
end

save HEt HEt;
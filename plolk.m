
xf=1:10:20000;
Nxf=length(xf);
LLk=zeros(1,Nxf);
for k=1:Nxf   
   s(3).time=xf(k);
   s=WorkVars(s,Root,LEAVES,NODES);
   TOPOLOGY=0;
   s=MarkRcurs(s,Root,[LEAVES,NODES],mu,TOPOLOGY);
   LLk(k)=LogLkd(s,Root,mu);
end

figure(201);
plot(xf,exp(LLk-max(LLk)));

f=pd(s(Root).ActI,mu,xf);
f=exp(f-max(f));
hold on; plot(xf,f,'--r');

nf=npd(s(Root).ActI,mu,xf);
nf=exp(nf-max(nf));
hold on; plot(xf,nf,':g');
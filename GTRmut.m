function [wordseq,langs,cogs,NS,L]=GTRmut(mu,lambda,s)

GlobalSwitches;

vocab=lambda/mu;
L=10*round(vocab);

x=vocab/L;
BF=[1 x];
BF=BF./sum(BF)

Q=[-1 1;1 -1]*diag(BF);
Q=Q-diag(diag(Q));
Q=Q-diag(sum(Q'));
Q=Q./abs(BF*diag(Q));

%align the model so that 2->1 at the same
%rate the Dollo-model would kill cognates

rate_21=mu/Q(2,1);

leaves=find([s.type]==LEAF);
NS=length(leaves);
N=2*NS-1;
Root=find([s.type]==ROOT);

seq=ones(N,L);
for site=1:L
   seq(Root,site)=disample(BF);
end
seq=GTRsim(s,Root,seq,L,rate_21,Q);

for k=1:N
   wordseq{k}=find(seq(k,:)==2);
end

langs = {s(leaves).Name}';
cogs=num2cell(strjust(num2str([1:L]'),'left'),2);

'hi I`m GTRmut'
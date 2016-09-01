function [seq,langs,cogs,NS,L,wdata]=MutationsRH(mu,p,lambda,s,bstd)
% function [seq,langs,cogs,NS,L,wdata]=MutationsRH(mu,p,lambda,s,bstd)
% mu is mean death rate
% p is cladagenic loss 
% lambda is word birth rate
% s is tree
% bstd is multiplier of std dev of word death rate so that
% total word death rate is ~LogNormal(mean=mu,std = mu*std)
GlobalSwitches;

leaves=find([s.type]==LEAF);
NS=length(leaves);
Root=find([s.type]==ROOT);

% see if we have rate heterogenity
if isempty(bstd) || bstd == 0
    mur=mu;
else
 %   mur=randg(v,v/mu);
  mur = randln(mu,mu*bstd,1);
%  disp(sprintf('At Root mu(edge)/mu(mean) = %g ',mur/mu));
end
n=MyPoisson(lambda/mur);

seq{Root}=1:n;
L=n;
wdata=[1:n;Root.*ones(1,n);zeros(1,n)];
[seq,L,wdata]=SeqSim(seq,L,wdata,s,Root,mu,p,lambda,bstd);

% langs=num2cell(strjust(num2str([leaves]'),'left'),2);
% DAVID 22 OCT 9:40
langs = {s(leaves).Name}';
cogs=num2cell(strjust(num2str([1:L]'),'left'),2);

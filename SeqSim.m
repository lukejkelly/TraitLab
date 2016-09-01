function [seq,eL,wdata]=SeqSim(seq,eL,wdata,s,k,mu,p,lambda,bstd)

%function [seq,eL]=SeqSim(seq,eL,s,Root,mu,p,lambda,bstd)
%
%Simulate word births at rate lambda
%and word deaths at rate mu down the
%two (or zero if k is leaf) edges below k

for kc=s(k).child
   
   dt=s(k).time-s(kc).time;
   
   %simulate new words along edge <k,kc>
   m=MyPoisson(lambda*dt);
   times=dt*rand(1,m);
   
   %simulate cladagenic loss of words from node <k>
   [oldwords,n]=BranchVocab(seq{k},p);
   
   %create new words
   newwords=(eL+1):(eL+m);
   eL=eL+m;
   
   %form the arrays of words and times for <k,kc> 
   words=[oldwords,newwords];
   dt=[dt.*ones(1,n),times];
   wdata=[wdata,[newwords;kc.*ones(1,m);times]];
   
   %now, which words make it to kc?
   % find local muation rate
   if isempty(bstd) || bstd == 0 
       % no rate het
       mur=mu;
   else
       % draw death rate ~LogNormal(mean = mu, std = mu*bstd)
       mur=randln(mu,mu*bstd,1);
  %     disp(sprintf('mu(%d)/mu(mean) = %g ',kc,mur/mu));
   end
   prob=exp(-mur*dt);
   rv=rand(1,n+m);  
   seq{kc}=words(rv<prob);
   
   %simulate down edges below kc
   [seq,eL,wdata]=SeqSim(seq,eL,wdata,s,kc,mu,p,lambda,bstd);
   
end

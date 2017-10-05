function [s toolarge]=WUprune(s,i,mu,p,kappa,cat)

global LOSTONES LEAF IN MIX MCMCCAT;
toolarge=0;
tlc1=0;
tlc2=0;


if s(i).type>LEAF

  c1=s(i).child(1);
  c2=s(i).child(2);

  if s(c1).mark
     [s tlc1]=WUprune(s,c1,mu,p,kappa,cat);
  end
  if s(c2).mark
     [s tlc2]=WUprune(s,c2,mu,p,kappa,cat);
  end

  dt1=s(i).time-s(c1).time;
  ef1=exp( -mu*dt1 );
  ncat1=cat(c1);
  pef1=p*ef1*(1-kappa)^ncat1; %probability of survival between i and c1, including catastrophes
  %kef1=(1-(1-kappa)^ncat1)/kappa;

  dt2=s(i).time-s(c2).time;
  ef2=exp( -mu*dt2 );
  ncat2=cat(c2);
  pef2=p*ef2*(1-kappa)^ncat2; 
  %kef2=(1-(1-kappa)^ncat2)/kappa;

  pef1u=1-s(c1).u*pef1;%probability of no descendant on side 1
  pef2u=1-s(c2).u*pef2;

  s(i).u = 1 - pef1u * pef2u ; %probability of having at least one descendant
  
  if MCMCCAT && kappa>0
      kef1=(1-(1-kappa)^ncat1)/kappa;
      kef2=(1-(1-kappa)^ncat2)/kappa;
  else
      kef1=1;
      kef2=1;
  end
      

  if LOSTONES

     %XXX next line goes with new Tu pruning
     s(i).v = s(c1).v*pef1*pef2u + s(c2).v*pef1u*pef2; %probability of having exactly one surviving descendant
     

         PartialLamSum1=(1-kappa)^ncat1*(1-ef1);
         PartialLamSum2=(1-kappa)^ncat2*(1-ef2);
         s(i).CatLamInt=s(c1).CatLamInt + s(c2).CatLamInt + s(c1).Tu*kef1+ s(c2).Tu*kef2;


     s(i).Tu = s(i).u - s(i).v; %probability of having at least two surviving descendants
     if s(i).Tu<0
         s(i).Tu=0;
         %error('Error: in WUprune, negative value of Tu was rounded to 0.');
         toolarge=1;
     end
     s(i).LamInt = s(c1).LamInt + s(c2).LamInt + s(c1).Tu*PartialLamSum1 + s(c2).Tu*PartialLamSum2;


  else %LOSTONES==0
         PartialLamSum1=(1-kappa)^ncat1*(1-ef1);
         PartialLamSum2=(1-kappa)^ncat2*(1-ef2);
         s(i).CatLamInt=s(c1).CatLamInt + s(c2).CatLamInt + s(c1).u*kef1 + s(c2).u*kef2;

         
        s(i).LamInt = s(c1).LamInt + s(c2).LamInt + s(c1).u*PartialLamSum1 + s(c2).u*PartialLamSum2;
        
  end

%   da1=s(i).ActI{1}; %deactivated on side 1: only descendants are on side 2
%   s(i).w(da1)=pef2*pef1u*s(c2).w(da1);
% 
%   da2=s(i).ActI{2};
%   s(i).w(da2)=pef1*pef2u*s(c1).w(da2);
% 
%   ab=s(i).ActI{3}; %descendants on both sides
%   s(i).w(ab)=pef1*pef2*(s(c1).w(ab).*s(c2).w(ab));
  
%   a:activated; da=deactivated; q=question marks
%   s(i).d is the number of descendant leaves

  s(i).d=s(c1).d+s(c2).d;
  
  % s(i).n is the probability that the trait was absent in all the
  % descendants of node i before missing data kicked in
  
  s(i).n=(pef1*s(c1).n+(1-pef1))*(pef2*s(c2).n+(1-pef2));
  
  da1a2=s(i).ActI{1};
  a1da2=s(i).ActI{2};
  a1a2=s(i).ActI{3};
  a1q2=s(i).ActI{4};
  q1a2=s(i).ActI{5};
  q1q2=s(i).ActI{6};
  q1da2=s(i).ActI{7};
  da1q2=s(i).ActI{8};
  
  s(i).w(a1da2)=pef1*s(c1).w(a1da2)*(1-pef2+pef2*s(c2).n);
  s(i).w(da1a2)=(1-pef1+pef1*s(c1).n).*pef2.*s(c2).w(da1a2);
  s(i).w(a1a2)=pef1*s(c1).w(a1a2).*pef2.*s(c2).w(a1a2);
  s(i).w(a1q2)=pef1*s(c1).w(a1q2).*(1-pef2 + pef2*s(c2).w(a1q2));
  s(i).w(q1a2)=(1-pef1 + pef1*s(c1).w(q1a2)).*pef2.*s(c2).w(q1a2);
  s(i).w(q1q2)=(1-pef1 + pef1*s(c1).w(q1q2)).*(1-pef2 + pef2*s(c2).w(q1q2));
  s(i).w(q1da2)=(1-pef1 + pef1*s(c1).w(q1da2))*(1-pef2+pef2*s(c2).n);
  s(i).w(da1q2)=(1-pef1+pef1*s(c1).n).*(1-pef2 + pef2*s(c2).w(da1q2));
  
  
  
  
  %PD is P[Ma=ma] for words born at regular birth events.
  s(i).PD(s(i).difCovI)=0;

  c1Cov=s(c1).CovI;


  PartialPD1=(1-kappa)^ncat1*(1-ef1);
  s(i).PD(c1Cov)=s(c1).PD(c1Cov)+s(c1).w(c1Cov).*PartialPD1;


  c2Cov=s(c2).CovI;


  PartialPD2=(1-kappa)^ncat2*(1-ef2);
  s(i).PD(c2Cov)=s(c2).PD(c2Cov)+s(c2).w(c2Cov).*PartialPD2;

  %PDcat is P[Ma=ma] for words born at catastrophic events.
  s(i).PDcat(s(i).difCovI)=0;
%   PartialPDcat1(c1Cov)=0;
%   PartialPDcat2(c2Cov)=0;

%  if DEPNU
      s(i).PDcat(c1Cov)=s(c1).PDcat(c1Cov) + kef1*ef1*s(c1).w(c1Cov);
      s(i).PDcat(c2Cov)=s(c2).PDcat(c2Cov) + kef2*ef2*s(c2).w(c2Cov);
%   else
%       for t=s(c1).cat
%          PartialPDcat1(c1Cov)=(1-kappa)*PartialPDcat1(c1Cov) + exp(-mu*(t-s(c1).time))*s(c1).w(c1Cov);
%       end
%       s(i).PDcat(c1Cov)=s(c1).PDcat(c1Cov) + PartialPDcat1(c1Cov);
% 
%       for t=s(c2).cat
%          PartialPDcat2(c2Cov)=(1-kappa)*PartialPDcat2(c2Cov) + exp(-mu*(t-s(c2).time))*s(c2).w(c2Cov);
%       end
%       s(i).PDcat(c2Cov)=s(c2).PDcat(c2Cov) + PartialPDcat2(c2Cov);
%   end



elseif s(i).type==LEAF

  if ~isempty(s(i).dat)
     s(i).w=((s(i).dat==IN)|(s(i).dat==MIX & s(i).xi<1)); %GKN Feb 08 added "& XI<1"
%     s(i).m=(s(i).dat==MIX);
  else
     s(i).w=[];
 %    s(i).m=[];
  end
  if LOSTONES
%      %XXX next line goes with new Tu pruning
%      if 0
%          s(i).Tw(:)=0;
%          s(i).Tw(i)=1;
%      end
     %XXX next line goes with new Tu pruning
     s(i).v=s(i).xi; 
  end
  s(i).u=s(i).xi;
  s(i).Tu=0;
  s(i).LamInt=0;
  s(i).CatLamInt=0;
  s(i).n=0;
  s(i).d=1;
else

  disp('unknown node type in WUprune.m'); keyboard;pause;

end

toolarge=toolarge || tlc1 || tlc2;

%if rand<10e-5, keyboard; end
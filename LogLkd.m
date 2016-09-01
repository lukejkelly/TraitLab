function [llkd,X]=LogLkd(state,lambda)

%call LogLkd(state) to get the integrated LL used in Markov() etc
%call LogLkd(state,[]) or LogLkd(state,lambda) to compute log(likelihood)
%assumes 'state' is properly set up - if really necessary, run
%state=MarkRcurs(state,state.nodes,true) to clean the state.
%
%At present it doesnt return the lambda value used in the former [] case.
%this might make sense (it would make the Lambda() function redundant)
%but would involve a trawl through all locations calling LogLkd() to correct the
%number of outputs. GKN 4 Jan 07


global LOSTONES MCMCCAT MISDAT;


s=state.tree;
Root=state.root;

Adam=s(Root).parent;
s(Adam).u=0;
if LOSTONES
  s(Adam).LamInt = s(Root).LamInt + s(Root).Tu;
else
  s(Adam).LamInt = s(Root).LamInt + s(Root).u;
end

%s(Adam).CatLamInt=s(Root).CatLamInt;
%s(Adam).w should stay equal to zeros(1,L)

ne=[s(Root).ActI{:}];
nd=length(ne);

s(Adam).PD(ne)=s(Root).PD(ne)+s(Root).w(ne);
%s(Adam).PDcat(ne)=s(Root).PDcat(ne);
%PDFromAdam(ne)=s(Root).w(ne)/(state.mu+state.rho*state.kappa);
%if MCMCCAT
    X=(s(Adam).LamInt + state.kappa*s(Root).CatLamInt)/state.mu;
%else
%    X=s(Adam).LamInt/state.mu;
%end

if MISDAT
    MissingNorm=sum([state.tree(state.leaves).nmis].*log(1-[state.tree(state.leaves).xi]))+sum((state.L-[state.tree(state.leaves).nmis]).*log([state.tree(state.leaves).xi]));
else
    MissingNorm=0;
end



if nargin==1 || isempty(lambda)
   %here is the calculation we usually make (for the MCMC and by default throughout)
   %we should not call the following the loglkd as we integrate lambda out with 1/lambda prior
   if nd>0
       %TODO there may be some erroneous factors of 2 here ?
       %if MCMCCAT
           llkd = MissingNorm+sum(log(s(Adam).PD(ne) + state.kappa*s(Root).PDcat(ne))) - nd*log(s(Adam).LamInt + state.kappa*s(Root).CatLamInt); %TODO check this against nargin==2 RJR 21/03/07
           %llkd = MissingNorm+sum(log(s(Adam).PD(ne) + state.kappa*s(Root).PDcat(ne))) - log(s(Adam).LamInt + state.kappa*s(Root).CatLamInt); %TODO check this against nargin==2 RJR 21/03/07
       %else
       %    llkd = MissingNorm+sum(log(s(Adam).PD(ne))) - nd*log(s(Adam).LamInt); %TODO check this against nargin==2 RJR 21/03/07
       %end
   else
       %disp('nd=0 in LogLkd.m');
       llkd = 0; %TODO 4/1/07 this is really a bug - nd=0 could mean, we looked and there were no cognates
       %or it could mean we didnt look. llkd=0 means we didnt look, but
       %that isnt how people expect LogLkd to behave. GKN
       %
       % Actually, this should not be fixed: nd=0 will never happen unless
       % we are trying to sample the prior, in which case llkd=0 is
       % correct. RJR 23/04/09
  end 

elseif nargin==2
   %this is the log likelihood, which we only use when we need the real thing for debugging
  % if nd==0
  %     llkd=0; %TODO 4/1/07 bug see comment above
       %disp('nd=0 in LogLkd.m');
  % else
        %llkd=-(lambda/state.mu)*s(Adam).LamInt - state.nu*s(Root).CatLamInt + sum(log(s(Adam).PD(ne)+state.kappa*s(Root).PDcat(ne))) + nd*log(lambda/state.mu) -gammaln(nd+1);
        llkd=MissingNorm-(lambda/state.mu)*s(Adam).LamInt - state.nu*s(Root).CatLamInt + sum(log(lambda/state.mu*s(Adam).PD(ne)+state.nu*s(Root).PDcat(ne))) ;%-gammaln(nd+1);
   %end
end

% if llkd>0
%     disp('Error in LogLkd: Log lkd is positive');
%     keyboard;
% end

%if llkd==-Inf
%    disp('Error in LogLkd: LogLkd=-Inf');keyboard;
%end

if imag(llkd) %|| isinf(llkd))
    disp('Error in LogLkd: Log lkd has imaginary part or is infinite');
    keyboard;
end


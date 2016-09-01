function [logprior,n]=LogPrior(prior,state)
global YULE FLAT ROOT

s=state.tree;
Root=state.root;
tl=state.length;

if prior.type==YULE
    %tl=TreeLength(s,Root);
    M=state.NS-1;
    logprior = gammaln(M+1) - M*log(tl) - log(M);
elseif prior.type==FLAT
    height=s(Root).time-min([s.time]);
    if height>prior.rootmax
        logprior=-inf;
    else
        %keyboard;
        if prior.isclade && ~isempty(prior.upboundclade)
            n=freeprogeny(s,[state.leaves,state.claderoot(prior.upboundclade)],Root);
            n(1)=[];
            if isempty(n)
                logprior=0;
            else
                logprior=-sum(log(s(Root).time-n)); %uniform
            end
            %logprior=-Nup*log(s(Root).time-mx); %uniform
        else
            logprior=-(state.NS-2)*log(height); %uniform
        end
    end
else
    disp('PRIOR type not recognised in  LogPrior');keyboard;pause;     
end

% %Taking into account catastrophies
% Ncat=state.ncat;
% %tl=TreeLength(s,Root);
% logprior=logprior - state.rho*tl +Ncat*log(state.rho*tl) -gammaln(Ncat+1);

for node=1:2*state.NS %[state.nodes state.leaves]
    if state.tree(node).type <ROOT
        dt=state.tree(state.tree(node).parent).time-state.tree(node).time;
        Ncat=state.cat(node);
        logprior=logprior -state.rho*dt +Ncat*log(state.rho*dt) -gammaln(Ncat+1);
    end
end


%Bug fix 14/12/04
%I used

%> int( theta^(2*n-3)*exp(-theta*X),theta=0..infinity);

%when I should have used

%> int( theta^(n-2)*exp(-theta*X),theta=0..infinity);
%                      1/(X^n)*X*GAMMA(n)/(n-1)

%so the matlab 

%M=2*state.NS-2;
%logprior = gammaln(M+2) - M*log(tl) - log( M*(M+1) );

%was wrong it should be

%M=state.NS-1;
%logprior = gammaln(M+1) - M*log(tl) - log(M);

%as above.

function [logprior,n]=LogPrior(prior,state)
global YULE FLAT ROOT TOPO LABHIST MCMCCAT LEAF BORROWING VARYRHO

s=state.tree;
Root=state.root;
tl=state.length;

if prior.type==YULE
    tl=TreeLength(s,Root);
    M=state.NS-1;
    logprior = gammaln(M+1) - M*log(tl) - log(M);
elseif prior.type==FLAT
    height=s(Root).time-min([s.time]);
    if height>prior.rootmax %| height<prior.rootmin
        logprior=-inf;
    else
        %GKN change 22/3/11 to add treatment of adamrange clades (to rootrange
        %clades)
        if prior.isclade && ~isempty(prior.upboundclade)
            tops=state.claderoot(prior.upboundclade);
            if any(prior.isupboundadamclade)
                ubac=prior.upboundclade(prior.isupboundadamclade); %upper bounded adam clades
                tops(prior.isupboundadamclade)=[s(state.claderoot(ubac)).parent]; %the bounded node is the parent
            end
            if any(tops==Root)
                logprior=0;
            else
                [n,v]=freeprogeny(s,[state.leaves,tops],tops,Root,prior);
                n(1)=[];
                if isempty(n)
                    logprior=0;
                else
                    %GKN change 23/3/11 to make uniform when there are
                    %adamrange nodes (added extra normalising term)
                    logprior=-sum(log(s(Root).time-n))+sum(log(prior.rootmax-n));
                end
            end
            %logprior=-Nup*log(s(Root).time-mx); %uniform
        else
            logprior=-(state.NS-2)*log(height); %uniform
        end
    end
else
    disp('PRIOR type not recognised in LogPrior');keyboard;pause;
end

if prior.topologyprior==TOPO
    % this is now correct - (was no '-1' but this error was masked by
    % an error in MCMC move 2, so was in fact the best value in
    % that setting) GKN 31/3/11
    logprior=logprior+sum(log([state.tree(state.nodes).d]-1));
    %elseif prior.topologyprior==LABHIST %% For now, assume that LABHIST is default
    %    nothing to do
    %else
    %    disp('Topology prior type not recognised in LogPrior');keyboard;pause;
end

% %Taking into account catastrophies
if MCMCCAT
    % Prior on number of catastrophes per branch is Poisson(rho)
    % If VARYRHO then we integrate rho analytically to get a Negative Binomial,
    % otherwise rho is fixed
    % If BORROWING then locations matter and these are uniform given counts
    if BORROWING
        if ~VARYRHO
            % Poisson counts and uniform locations
            % Prior on number of catastrophes k on and positions w on each
            % branch is k ~ Poisson(rho * dt) and w | k ~ k! / (dt)^k so
            % mass/density function is pi(w, k | rho) = rho^k exp(-rho * dt)
            for node = 1:(2 * state.NS)
                if state.tree(node).type < ROOT
                    dt       = state.tree(state.tree(node).parent).time - state.tree(node).time;
                    Ncat     = state.cat(node);
                    logprior = logprior - state.rho * dt + Ncat * log(state.rho);
                end
            end
        else
            % Negative Binomial counts and uniform locations on entire tree
            % a and b from LogRhoPrior
            a = 1.5; b = 5e3;
            logprior = logprior + gammaln(a + state.ncat) - (a + state.ncat) * log(state.length + b);
        end
    else
        % We only consider number of catastrophes on each branch, not locations
        if ~VARYRHO
            % Poisson prior on branch counts
            for node = 1:2*state.NS
                if state.tree(node).type <ROOT
                    dt       = state.tree(state.tree(node).parent).time - state.tree(node).time;
                    Ncat     = state.cat(node);
                    logprior = logprior - state.rho * dt + Ncat * log(state.rho * dt) - gammaln(Ncat + 1);
                end
            end
        else
            % Negative Binomial counts, a and b from LogRhoPrior
            a = 1.5; b = 5e3;
            logprior = logprior + gammaln(a + state.ncat) - (a + state.ncat) * log(state.length + b);
            for node = 1:2*state.NS
                if state.tree(node).type <ROOT
                    dt       = state.tree(state.tree(node).parent).time - state.tree(node).time;
                    Ncat     = state.cat(node);
                    logprior = logprior + Ncat * log(dt) - gammaln(Ncat + 1);
                end
            end
        end
    end
end

%Bug fix september 07

%I didnt implement what I wrote in the paper, and what I wrote was better
%so that's now here - I was taking the difference between the root and the
%lowest node free node, now it is the root and the relevant lower bound.

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

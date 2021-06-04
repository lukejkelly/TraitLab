function state=makestate(prior,mu,lambda,p,rho,kappa,content,s,ListOfCats, beta)

global MIX DONTMOVECATS ANST BORROWING

% To avoid conflicts with non-borrowing code.
if nargin <= 9, beta = 0; end

%function state=makestate(prior,mu,p,content,s)
%
%typically used to build start state for MCMC

%%
%state, data
state=tree2state(s);

%TODO check that this doesn't lead to any problems; should ListOfCats
%always be something? RJR 28/02/07
%GKN added isempty Jan 08
% LJK added all(ListOfCats == 0) for true state when taxa dropped.
if nargin<=8 || isempty(ListOfCats) || all(ListOfCats == 0)
    ListOfCats=zeros(2*state.NS,1);
    if DONTMOVECATS
        for i=1:2*state.NS
            if s(i).type<=ANST, ListOfCats(i)=1; end
        end
    end
end

if state.NS~=content.NS
   disp('unequal # taxa found in InstallData args content and state.tree');keyboard;pause;
end


%write data onto tree leaves
[state.tree, ok] = mergetreedata(state.tree,content.array,content.language);
if ~ok
    dbstack,disp('Problem writing data into tree in makestate');keyboard;pause; %GKN 27/10
end

%%
%clades
state.claderoot=[];
if prior.isclade
   state.tree = treeclades(state.tree,prior.clade);
   state=UpdateClades(state,[state.leaves,state.nodes],size(prior.clade,2));
end

%%
%state parameters and properties of state
state.mu=mu;
state.p=p;
state.lambda=lambda;
state.rho=rho;
state.kappa=kappa;
state.nu=state.lambda*state.kappa/state.mu;
state.L=content.L;
state.cat=ListOfCats(:); % LUKE 05/09/2016 to ensure this appears as a column vector.
state.ncat=sum(state.cat);
state.length=TreeLength(state.tree,state.root);
state.beta = beta; % LUKE

%%
%initial value of xi
for i=state.leaves
    state.tree(i).nmis=length(find(state.tree(i).dat==MIX));
    state.tree(i).xi=1-state.tree(i).nmis/state.L;
    if state.tree(i).xi==1, state.tree(i).xi=.9999; end
end

%%
%5) propagate the data up the tree, computing work array values
state.tree=WorkVars(state.NS,state.L,state.tree,state.root);
TOPOLOGY=1;
state = MarkRcurs(state,[state.leaves,state.nodes],TOPOLOGY);

%%
%maintain the tree log-likelihood and log-prior
% state.loglkd=LogLkd(state);

% Setting initial catastrophe locations in borrowing model when we start from an
% old state.
if BORROWING
    % Catastrophes cannot go on the edge <Adam, Root>.
    for i = (find(state.cat & ([state.tree.type]' <= ANST) ...
            & cellfun(@isempty, {state.tree.catloc})'))'
        % Catastrophes are uniformly distributed along edge <pa(i), i>.
        % LJK 04/06/21 sorting
        state.tree(i).catloc = sort(rand(state.cat(i), 1));
    end
end

if BORROWING    % Luke 28/01/2014
    [state.loglkd, state.fullloglkd] = logLkd2(state);
else
    state.loglkd = LogLkd(state);
    state.fullloglkd = LogLkd(state,state.lambda);
end

state.logprior=LogPrior(prior,state);

%assign lambda to a random draw from the posterior with 1/lambda prior
%state.lambda=Lambda(state);

%%
% update fullloglkd (i.e., without lambda integrated out)
% state.fullloglkd = LogLkd(state,state.lambda); % Luke 28/01/2014

%%
% Clear persistent variables in SDLT code. LUKE 04/09/2016
clear logLkd2_m patternCounts patternMeans

function [nstate, U, OK, logq]=DelCat(state)
% [nstate, U, OK, logq]=DelCat(state)
% RJR, last modified on 09/05/07
% LJK, last modified on 07/02/14

global BORROWING

nstate=state;
OK=1;
if state.ncat==0
    OK=0;
    U=[];
else
    r=ceil(state.ncat*rand);
    node=find(cumsum(state.cat)>=r,1);
    nstate.cat(node)=state.cat(node)-1;
    nstate.ncat=state.ncat-1;

    % Luke 07/02/14
    % Picking one of the catastrophes on <node, pa(node)> and getting rid of it.
    if BORROWING
    	nstate.tree(node).catloc( ceil(length(nstate.tree(node).catloc) * rand) ) = [];
    end

    if state.cat(node)~=0
            %t=ceil(numel(nstate.tree(node).cat)*rand);
            %nstate.tree(node).cat=[state.tree(node).cat(1:(t-1)), state.tree(node).cat((t+1):end)];
            U=above(node, state.tree, state.root);
    else
        OK=0;
        U=[];
        disp('Error in DelCat: Attempted to delete a catastrophe from a branch without catastrophes.');
        keyboard;
    end   
end

if OK
    dt=nstate.tree(nstate.tree(node).parent).time-nstate.tree(node).time;
    %logq=  +log(state.ncat)   -log(dt) +log(2*nstate.NS-2);
    % logq=-log(state.length) +log(dt) +log(state.ncat) -log(state.cat(node)); Luke 19/02/14 put this inside below if-statement.

    % Luke 19/02/14 - added this if-statement to account for the fact that proposal distribution with borrowing applies to
    % catastrophes at specific times. Catastrophes on a node can be identified by their times. Proposal distribution to remove
    % a specific catastrophe  on <pa(j), j> is 1 / nstate.ncat. Proposal distribution to add one again in the same location
    % is 1 / state.length. 
    if BORROWING
        logq = log(state.ncat) - log(state.length);
    else
        logq=-log(state.length) +log(dt) +log(state.ncat) -log(state.cat(node));
    end

else
    logq=0;
end

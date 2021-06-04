function [nstate, U, OK, logq]=AddCat(state)
% [state, U, OK, logq]=AddCat(state)
% RJR. Last modified on 07/05/07

global ADAM ROOT BORROWING

nstate=state;
L=2*state.NS; %number of nodes
OK=0;
roottime=state.tree(state.root).time;

while ~OK
    r=ceil(L*rand);
    if any(state.tree(r).type==[ADAM ROOT])
        OK=0;
    else
        dt=state.tree(state.tree(r).parent).time - state.tree(r).time;
        if rand<dt/roottime
            OK=1;
        else
            OK=0;
        end
    end
end


%dt=state.tree(state.tree(r).parent).time - state.tree(r).time;
% t=state.tree(r).time + rand*dt;
% nstate.tree(r).cat=InsertOrdered(state.tree(r).cat, t);
nstate.cat(r)=state.cat(r)+1;
nstate.ncat=state.ncat+1;

% Luke 07/02/2014 - adding a U[0, 1] rv to <r, pa(r)>; 04/06/2021 - sorting
if BORROWING
    nstate.tree(r).catloc = sort([ nstate.tree(r).catloc, rand ]);
end

U=above(r, nstate.tree, nstate.root);
%logq=  -log(nstate.ncat) +log(dt) +log(L-2);
% logq=log(state.length)-log(dt)-log(nstate.ncat)+log(nstate.cat(r)); % Luke 19/02/14 put this inside below if-statement.

% Luke 19/02/14 - added this if-statement to account for the fact that proposal distribution with borrowing applies to
% catastrophes at specific times. Catastrophes on a node can be identified by their times. Proposal distribution to
% add a catastrophe at time t on <j, pa(j)> is 1 / state.length. Proposal distribution to remove the catastrophe again
% is simply 1 / nstate.ncat.
if BORROWING
    logq = log(state.length) - log(nstate.ncat);
else
    logq=log(state.length)-log(dt)-log(nstate.ncat)+log(nstate.cat(r));
end

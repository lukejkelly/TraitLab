function [nstate,U,OK,logq]=MoveCat(state)
% [nstate,U,OK,logq]=MoveCat(state)
% Moves a catastrophe (chosen UAR) to a legal neighbouring node: child, parent or sibling
% RJR, 07/05/07
% LJK 07/02/14, account for catastrophes in the presence of borrowing
% LJK 04/06/21, sorting catastrophe locations
% LJK 19/07/21, new catastrophe location on branch

global BORROWING

nstate=state;

if state.ncat==0
    OK=0;
    U=[];
    logq=0;
else
    OK=1;

    r=ceil(rand*state.ncat);
    old=find(cumsum(state.cat)>=r,1);

    [new, q1, q2]=GetLegal(state.tree,old,state.root);

    nstate.cat(old)=state.cat(old)-1;
    nstate.cat(new)=state.cat(new)+1;

    if BORROWING
        % Moving catastrophe location
    	ind = ceil(length(nstate.tree(old).catloc) * rand);
    	nstate.tree(old).catloc(ind) = [];
        nstate.tree(new).catloc = sort([nstate.tree(new).catloc, rand]);
    end

    % Moving a catastrophe is equivalent to an add/delete with the condition
    % that we land on a nearby branch. We choose a catastrophe with probability
    % 1 / N. The catastrophe is on branch i and we then choose a neighbour j
    % with probability 1 / q_i. Finally we choose a location with density
    % 1 / dt_j. The 1 / N terms cancel.
    if BORROWING
        dt1 = nstate.tree(nstate.tree(old).parent).time - nstate.tree(old).time;
        dt2 = nstate.tree(nstate.tree(new).parent).time - nstate.tree(new).time;
        logq = log(q1) - log(q2) - log(dt1) + log(dt2);
    else
        logq=+log(q1)-log(q2)-log(state.cat(old))+log(nstate.cat(new));
    end

    U=above([old,new],nstate.tree,nstate.root);
end

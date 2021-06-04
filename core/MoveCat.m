function [nstate,U,OK,logq]=MoveCat(state)
% [nstate,U,OK,logq]=MoveCat(state)
% Moves a catastrophy (chosen UAR) to a legal neighbouring node: child, parent or sibling
% RJR, 07/05/07

% LJK, 07/02/14 - to account for catastrophes in the presence of borrowing.

global BORROWING % Luke 19/02/2014

nstate=state;

if state.ncat==0
    OK=0;
    U=[];
    logq=0;
else
    OK=1;

    r=ceil(rand*state.ncat);
    old=find(cumsum(state.cat)>=r,1);

%    if length(old)==0, keyboard; end

    [new, q1, q2]=GetLegal(state.tree,old,state.root);

    del=ceil(rand*numel(state.cat(old)));
    %nstate.tree(old).cat=[nstate.tree(old).cat(1:del-1) nstate.tree(old).cat(del+1:end)];
    nstate.cat(old)=state.cat(old)-1;
    nstate.cat(new)=state.cat(new)+1;
    %dt1=nstate.tree(nstate.tree(old).parent).time - nstate.tree(old).time;
    dt2=nstate.tree(nstate.tree(new).parent).time - nstate.tree(new).time;
    %nstate.tree(new).cat=InsertOrdered(nstate.tree(new).cat,rand*dt2+nstate.tree(new).time);

    % Luke 07/02/14
    % Picking one of the catastrophes on <old, pa(old)> and adding it to <new, pa(new)>.
    if BORROWING

        % Step one - picking which catastrophe to move.
    	ind = ceil(length(nstate.tree(old).catloc) * rand);

        % Step two - placing the catastrophe on the 'new' branch and removing it from the 'old' one.
        % LJK 04/06/21 sorting
    	nstate.tree(new).catloc = sort([ nstate.tree(new).catloc, ...
                                         nstate.tree(old).catloc(ind) ]);
    	nstate.tree(old).catloc(ind) = [];

    end

    %logq=log(dt2)-log(dt1)+log(q1)-log(q2);
    % logq=+log(q1)-log(q2)-log(state.cat(old))+log(nstate.cat(new)); Luke 19/02/14 put this inside below if-statement.

    % Luke 11/5/16
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

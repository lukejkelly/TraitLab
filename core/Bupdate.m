function [state,U,TOPOLOGY]=Bupdate(state,i,j,k,newage)

% [state,U,TOPOLOGY]=Bupdate(state,i,j,k,newage)
% Reconnects an edge into another edge.
% k is j's parent; iP is i's parent
% iP becomes k's child and j's parent, with time newage

global OTHER ANST ROOT BORROWING % Luke 09/02/14 --- catastrophes and borrowing.

s=state.tree;
Root=state.root;

iP=s(i).parent;
PiP=s(iP).parent;
CiP=s(iP).child(OTHER(s(i).sibling));

s(PiP).child(s(iP).sibling)=CiP;
s(CiP).parent=PiP;
s(CiP).sibling=s(iP).sibling;

s(j).parent=iP;
s(k).child(s(j).sibling)=iP;
s(iP).sibling=s(j).sibling;
s(iP).child(OTHER(s(i).sibling))=j;
s(j).sibling=OTHER(s(i).sibling);
s(iP).parent=k;

% Luke --- If iP is the root then after moving <iP, i>, j is the new root, so we
% remove the catastrophes from CiP and place them on iP. The total number of
% catastrophes is preserved. We perform the reverse move if j is the root.

if BORROWING

    if iP == Root

        % iP becomes an internal node and CiP the new root.
        s(Root).type = ANST;
        Root = CiP;
        s(CiP).type = ROOT;

        % If iP is the root then it has no catastrophes. We move the catastrophes on
        % CiP to iP as CiP is the new root.
        state.cat(iP) = state.cat(CiP);
        s(iP).catloc = s(CiP).catloc;

        state.cat(CiP) = 0;
        s(CiP).catloc = [];

    elseif j == Root

        % j becomes an internal node and iP the new root.
        s(Root).type = ANST;
        Root = iP;
        s(iP).type = ROOT;

        % If j is the root then iP becomes the new root. We move the catastrophes on
        % iP to j and set the number of catastrophes on iP to be zero.
        state.cat(j) = state.cat(iP);
        s(j).catloc = s(iP).catloc;

        state.cat(iP) = 0;
        s(iP).catloc = [];

    end

    % If neither iP nor j is the root, then the number of catastrophes on each branch
    % remains constant.

else % How catastrophes are treated when there is no borrowing.

    if iP==Root
        s(Root).type=ANST;
        Root=CiP;
        s(Root).type=ROOT;
        state.ncat=state.ncat-state.cat(Root);
        state.cat(Root)=0;
    elseif j==Root
        s(Root).type=ANST;
        Root=iP;
        s(Root).type=ROOT;
        state.ncat=state.ncat-state.cat(Root);
        state.cat(Root)=0;
    end

    % The following two lines have been moved outside the if/else statement.
    % oldage=s(iP).time;
    % s(iP).time=newage;

    % Update the catastrophes
    state.cat(CiP)=state.cat(CiP)+state.cat(iP);
    state.cat(iP)=0;

    % ti=s(i).time;
    % s(i).cat=(s(i).cat-ti)/(oldage-ti)*(newage-ti)+ti;
    % 
    % tc=find(s(j).cat>newage);
    % if numel(tc)>0
    %     s(iP).cat=s(j).cat(tc);
    %     state.cat(find(state.cat==j,numel(tc)))=iP;
    %     %s(j).cat=s(j).cat(find(s(j).cat<=newage));
    %     s(j).cat(tc)=[];
    % end

    state.cat(iP)=binornd(state.cat(j),(newage-s(j).time)/(s(k).time-s(j).time));
    state.cat(j)=state.cat(j)-state.cat(iP);

end

oldage = s(iP).time;
s(iP).time = newage;

U=above([PiP,iP],s,Root);

state.tree=s;
state.root=Root;

if any([iP, PiP, CiP, j]==Root) %just to make sure
    state.length=TreeLength(state.tree,Root);
else
    state.length=state.length-oldage+newage;
end
TOPOLOGY=1;


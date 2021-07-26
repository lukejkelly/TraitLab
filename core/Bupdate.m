function [state, U, TOPOLOGY] = Bupdate(state, i, j, k, newage, ncat, cat, loc)

% [state, U, TOPOLOGY]=Bupdate(state, i, j, k, newage)
% Reconnects an edge into another edge
% k is j's parent; iP is i's parent
% iP becomes k's child and j's parent, with time newage

% LJK 09/02/14: catastrophes and borrowing
% LJK 21/07/21: same catastrophe definition, account for resampling

global OTHER ANST ROOT MCMCCAT

% Update catastrophe count in advance of SPR move
if MCMCCAT
    state = Bcats.updateCats(state, i, j, ncat, cat, loc);
end

s = state.tree;
Root = state.root;

iP = s(i).parent;
PiP = s(iP).parent;
CiP = s(iP).child(OTHER(s(i).sibling));

s(PiP).child(s(iP).sibling)=CiP;
s(CiP).parent = PiP;
s(CiP).sibling = s(iP).sibling;

s(j).parent = iP;
s(k).child(s(j).sibling)=iP;
s(iP).sibling = s(j).sibling;
s(iP).child(OTHER(s(i).sibling))=j;
s(j).sibling = OTHER(s(i).sibling);
s(iP).parent = k;

if iP == Root
    % iP becomes an internal node and CiP the new root
    s(Root).type = ANST;
    Root = CiP;
    s(Root).type = ROOT;
elseif j == Root
    % iP is new root
    s(Root).type = ANST;
    Root = iP;
    s(Root).type = ROOT;
end

oldage = s(iP).time;
s(iP).time = newage;

U = above([PiP, iP], s, Root);

state.tree = s;
state.root = Root;

if any([iP, PiP, CiP, j] == Root)
    state.length = TreeLength(state.tree, Root);
else
    state.length = state.length - oldage + newage;
end
TOPOLOGY = 1;

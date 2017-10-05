function [nstate,U,TOPOLOGY]=Supdate(state,i,newage)
% [nstate,U,TOPOLOGY]=Supdate(state,i,newage)
% mode node i to time newage
% GKN; last modified by RJR on 09/05/07


global LEAF ROOT ANST

nstate=state;
oldage=state.tree(i).time;
nstate.tree(i).time=newage;

% if state.tree(i).type<ROOT % i is LEAF or ANST, there are catastrophes above it
%     t2=state.tree(state.tree(i).parent).time;
%     nstate.tree(i).cat=(state.tree(i).cat-oldage)*(t2-newage)/(t2-oldage)+newage; %linear transformation
% end
% 
% if state.tree(i).type>LEAF % i is ANST or ROOT, there are catastrophes below it
%     c1=state.tree(i).child(1);
%     c2=state.tree(i).child(2);
%     tc1=state.tree(c1).time;
%     tc2=state.tree(c2).time;
%     nstate.tree(c1).cat=(state.tree(c1).cat-tc1)*(newage-tc1)/(oldage-tc1)+tc1;
%     nstate.tree(c2).cat=(state.tree(c2).cat-tc2)*(newage-tc2)/(oldage-tc2)+tc2;
% end  

switch state.tree(i).type
    case LEAF
        disp('Warning: Supdate should not be used to modify age of a leaf.');
        keyboard; pause;
        nstate.length=state.length+oldage-newage;
    case ANST
        nstate.length=state.length-oldage+newage;
    case ROOT
        nstate.length=state.length-2*oldage +2*newage;
end


U=above(i,state.tree,state.root);
TOPOLOGY=0;

function str=wnexcattree(s,node,cat)
%str=wnexcattree(s,node,cat)

global LEAF ANST ROOT ADAM;


% check what type of node we have
switch s(node).type
case ROOT
    str = ['(' wnexcattree(s,s(node).child(1),cat) ',' wnexcattree(s,s(node).child(2),cat) ')']; 
case ANST
    str = ['(' wnexcattree(s,s(node).child(1),cat) ',' wnexcattree(s,s(node).child(2),cat) ')' ':' num2str(cat(node)+.1)]; 
case LEAF
    str = strcat(s(node).Name,':',num2str(cat(node)+.1)); 
case ADAM
    disp('wnexcattree unexpectedly passed ADAM node')
    str = wnexcattree(s,s(node).child,cat);
otherwise 
    disp(['Error in wnexcattree - s.type not recognised at node ' num2str(node)])
    keyboard
end
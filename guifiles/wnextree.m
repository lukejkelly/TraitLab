function str = wnextree(s,node)

global LEAF ANST ROOT ADAM


% check what type of node we have
switch s(node).type
case ROOT
    str = ['(' wnextree(s,s(node).child(1)) ',' wnextree(s,s(node).child(2)) ')']; 
case ANST
    time = abs(s(s(node).parent).time - s(node).time);
    str = ['(' wnextree(s,s(node).child(1)) ',' wnextree(s,s(node).child(2)) ')' ':' num2str(time)]; 
case LEAF
    time = abs(s(s(node).parent).time - s(node).time);
    str = strcat(s(node).Name,':',num2str(time)); 
case ADAM
    disp('wnextree unexpectedly passed ADAM node')
    str = wnextree(s,s(node).child);
otherwise 
    disp(['Error in wnextree - s.type not recognised at node ' num2str(node)])
    keyboard
end
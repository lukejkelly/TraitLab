
clear;
disp('Full tree');
load run2;
draw(state.tree,1,0,'IEtree');
disp('Show cognates');
pause;
for k=1:10, draw(state.tree,1,0,'IEtree',k,''); pause; end
draw(state.tree,1,0,'IEtree');

disp('Next step prepare MCMC');
pause;
tree;

%clear;
%load run3;
%draw(state.tree,2,0,'IEtree 20');
 
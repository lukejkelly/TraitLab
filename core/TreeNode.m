function node=TreeNode(par,sib,chd,t,name,type)

%GlobalSwitches;

node=pop('node');

node.parent=par;
node.sibling=sib;
node.child=chd;
node.time=t;
node.Name=name;
node.type=type;

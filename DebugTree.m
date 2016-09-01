function sTR=DebugTree(NS,ThetaTR)

global ROOT LEAF ADAM

if NS~=2, disp('wrong # leaves for DebugTree');keyboard;pause; end

sTR(3)=TreeNode(4,1,[1,2],3000,num2str(3),ROOT);
sTR(2)=TreeNode(3,2,[],0,'lang2',LEAF);
sTR(1)=TreeNode(3,1,[],0,'lang1',LEAF);

sTR(4)=TreeNode([],[],3,realmax,'Adam',ADAM);
sTR(3).parent=4;
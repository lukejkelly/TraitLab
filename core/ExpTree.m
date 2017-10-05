function sTR=ExpTree(NS,ThetaTR,tms)
%What is tms? It isn't used. RJR 28/02/07

global ROOT ANST LEAF ADAM MISDAT

t=0;
n=1;
s(1)=TreeNode([],1,[],t(n),num2str(n),ROOT);

nl=2;
t=-(1/ThetaTR)*log(rand(1,2));
n=[2,3];
m=4;
p=[1,1];

while nl<NS
   
   [ct,I]=min(t);
   cn=n(I);
   s(cn)=TreeNode(p(I),[],[],t(I),num2str(n(I)),ANST);
             
   n(I)=[];
   p(I)=[];
   t(I)=[];
   
   n=[n,m,m+1];
   p=[p,cn,cn];      
   t=[t,ct-(1/ThetaTR)*log(rand(1,2))];   

   
   m=m+2;   
   nl=length(n);
   
end

[ct,I]=min(t);
for k=1:length(n)%this looks like a bug to me: with L
   s(n(k))=TreeNode(p(k),[],[],ct,num2str(n(k)),LEAF);
end

LEAVEStr=find([s.type]==LEAF);
NODEStr=find([s.type]>LEAF);
RootTR=1;

v=[max([s.time])-[s.time]];
N=2*NS-1;
for j=1:N
   s(j).time=v(j);
end
for k=1:N
   if s(k).type<ROOT
      kp=s(k).parent;
      cl=s(kp).child;
      s(kp).child=[cl,k];
      s(k).sibling=length(s(kp).child);
   else
      ;
   end
end

%relabel tree so NODES 1 to NS (incl) are LEAF and NS+1 to 2NS-1 are ANST 
for k=LEAVEStr
   MAP(k)=find(LEAVEStr==k);
end
for k=NODEStr
   MAP(k)=NS+find(NODEStr==k);
end
for k=[LEAVEStr,NODEStr]
   s(k).child=MAP(s(k).child);
   s(k).parent=MAP(s(k).parent);
   s(k).Name=num2str(MAP(k));
   sTR(MAP(k))=s(k);
end
RootTR=MAP(RootTR);

sTR(2*NS)=TreeNode([],[],RootTR,realmax,'Adam',ADAM);

sTR(RootTR).parent=2*NS;

% if MISDAT
%     %for k=1:NS
%     %    s(k).xi=1-rand^3;
%     %end
%     disp('Vector of missing data parameter values (xi):');
%     for k=LEAVEStr
%         s(k).Name
%         s(k).xi
%     end
% end

% give slightly meaningful names to leaves
% leafnames = sprintf('lang%5.0f',1:NS);
% leafnames = reshape(leafnames,9,NS)';
% leafnames = cellstr(leafnames);
% leafnames = strrep(leafnames,' ','');
% [sTR(1:NS).Name] = deal(leafnames{:});

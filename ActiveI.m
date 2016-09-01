function s=ActiveI(s,i)

%  Modified by RJR on 14/11/07 for missing data

global LOSTONES LEAF IN MIX
global za1 za2 zb1 zb2 zq1 zq2

% 1    01
% 2    10
% 3    11
% 4    1?
% 5    ?1
% 6    ??
% 7    ?0
% 8    0?

if s(i).type>LEAF
   
   c1=s(i).child(1);
   c2=s(i).child(2);

   if s(c1).mark==1
      s=ActiveI(s,c1);
   end
   if s(c2).mark==1
      s=ActiveI(s,c2);
   end
   
   a1=[s(c1).ActI{1:5}];
   a2=[s(c2).ActI{1:5}];
   q1=[s(c1).ActI{6:8}];
   q2=[s(c2).ActI{6:8}];
%   b1=[a1 q1];
%   b2=[a2 q2];
   
   %s(i).ActI{1}=setdiff(a2,b1); %in i from 2 not in 1
   %s(i).ActI{2}=setdiff(a1,b2); %in i from 1 not in 2
   %s(i).ActI{3}=intersect(a1,a2); %from both 1 and 2
   %s(i).ActI{4}=intersect(a1,q2)
   %s(i).ActI{5}=intersect(q1,a2)
   %s(i).ActI{6}=intersect(q1,q2)
   %s(i).ActI{7}=setdiff(q1,b2)
   %s(i).ActI{8}=setdiff(q2,b1)
   %the above is easier to read...
   
   %but the following is faster
%    za(:)=0;
%    zb(:)=0;
%    za(a1)=1;
%    zb(a2)=1;
   
   za1(:)=0;
   za2(:)=0;
   zq1(:)=0;
   zq2(:)=0;
   zb1(:)=0;
   zb2(:)=0;
   
   za1(a1)=1;
   za2(a2)=1;
   zq1(q1)=1;
   zq2(q2)=1;

   zb1(:)=~(za1 | zq1);
   zb2(:)=~(za2 | zq2);

   
   s(i).ActI{1}=find(za2&zb1);
   s(i).ActI{2}=find(za1&zb2);
   s(i).ActI{3}=find(za1&za2);
   s(i).ActI{4}=find(za1&zq2);
   s(i).ActI{5}=find(zq1&za2);
   s(i).ActI{6}=find(zq1&zq2);
   s(i).ActI{7}=find(zq1&zb2);
   s(i).ActI{8}=find(zq2&zb1);
   
   %%GKN 28/10 start
   %s(i).clade=union(s(c1).clade,s(c2).clade);
   %%GKN 28/10 end
      
   %TCovI cognates are present at just one leaf
   %so they are covered by a node if they are active there
   if LOSTONES
      s(i).TCovI=[s(c1).TCovI,s(c2).TCovI];
   end
   
   
elseif s(i).type==LEAF

   s(i).ActI{1}=[];
   s(i).ActI{2}=[];
   s(i).ActI{4}=[];
   s(i).ActI{5}=[];
   s(i).ActI{7}=[];
   s(i).ActI{8}=[];
   
   if ~isempty(s(i).dat)
      s(i).ActI{3}=find(s(i).dat==IN);
      s(i).ActI{6}=find(s(i).dat==MIX);
   else
      s(i).ActI{3}=[];
      s(i).ActI{6}=[];
   end
   
   %TCovI cognates are present at just one leaf.
   %They have a separate labeling, by their leaf #
   if LOSTONES
      s(i).TCovI=i;
   end
   
else

   disp('unknown node type in Active.m');keyboard;pause;
   
end


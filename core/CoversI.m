function s=CoversI(s,i)

%function [s,c]=CoversI(s,i,c)

global za1 zb1 zb2 LEAF

c1=s(i).child(1);
c2=s(i).child(2);

%s(c1).CovI=intersect(s(i).ActI{2},s(i).CovI);
%s(c2).CovI=intersect(s(i).ActI{1},s(i).CovI);
%childcover=union(s(c1).CovI,s(c2).CovI);
%s(i).difCovI=setdiff(s(i).CovI,childcover);
%the above is easier to read

%but the following is faster
za1(:)=0;
zb1(:)=0;
zb2(:)=0;

za1(s(i).CovI)=1;

zb1(s(i).ActI{2})=1;
zb1(s(i).ActI{4})=1;
s(c1).CovI=find(za1&zb1);

zb2(s(i).ActI{1})=1;
zb2(s(i).ActI{5})=1;
s(c2).CovI=find(za1&zb2);

s(i).difCovI=find(za1&~(zb1|zb2));





if s(c1).mark && ~isempty(s(c1).CovI) && s(c1).type>LEAF
   s=CoversI(s,c1);
   %[s,c]=CoversI(s,c1,c);
end

if s(c2).mark && ~isempty(s(c2).CovI) && s(c2).type>LEAF
   s=CoversI(s,c2);
   %[s,c]=CoversI(s,c2,c);
end

%%GKN 28/10 maintain a list (state.claderoot) of clade-mrca's 
%c(intersect(s(c1).clade,s(c2).clade))=i;

function X=above(U,s,Root)

% function X=above(U,s,Root)
%
% U: a list of node labels on s
% X: nodes on s above or equal U
% U is a subset of X, except that Adam is not in X even if it is in U
% CHANGE GKN 18/10/02

V=U;
W=[];
for k=1:length(U)
   p=U(k);
   if isempty(s(p).parent)
      V(k)=[];
   else
      while p~=Root
         p=s(p).parent;
         if s(p).mark==1
            break;
         else
            W=[W,p];
            s(p).mark=1;
         end
      end
   end
end

X=[V,W];

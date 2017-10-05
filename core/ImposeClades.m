function s=ImposeClades(s,r,c,filename)

% function s=ImposeClades(s,r,c)
%
% s: state.tree, r: state.claderoot, c: prior.clade
% Used in inisialisation of MCMC to make a tree
% satisfying clade age constraints
%
% assumes tree s actually has the clades in c
% otherwise this algorithm will fail
% TODO fix this - build clade topology 
% and leaf times into exptree

GlobalSwitches;

n=length(r);

for k=1:n
   
   rr=c{k}.rootrange;
   if ~isempty(rr)
      i=r(k);
      IsInRange=InRange(s(i).time,rr);
      if ~IsInRange
         if rr(1)==-inf
            target=rr(2);
         elseif rr(2)==inf
            target=rr(1);
         else
            target=mean(rr);
         end
         scale=target/s(i).time;
         b=[i,below(s,i)];
         for j=b
            if s(j).type~=LEAF
               s(j).time=s(j).time*scale;
            end
         end
      end
   end
   
   ar=c{k}.adamrange;
   if ~isempty(ar)
      i=s(r(k)).parent;
      IsInRange=InRange(s(i).time,ar);
      if ~IsInRange
         if ar(1)==-inf
            s(i).time=ar(2);
         elseif ar(2)==inf
            s(i).time=ar(1);
         else
            s(i).time=mean(ar);
         end
      end
   end
 
end

fn=[];
FAIL=0;
for k=1:size(s,2)
   for ck=s(k).child
      if s(k).time<s(ck).time
         FAIL=1;
         fn=[fn,[k;ck]]
      end
   end
end

if FAIL
   disp('ImposeClades.m failed on the following edges');
   disp(fn);
   keyboard;pause;
elseif nargin==4
   [str,ok] = stype2nexus(s,'Suitable MCMC init tree respecting clade age constraints','TREE',[]);
   fid = fopen(filename,'w');fprintf(fid,str);fclose(fid);
end

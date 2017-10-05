function s=FinishBuilding(s,i,B)

%FinishBuilding put parent info into a tree with just child info
%
%s=FinishBuilding(s,i)
%
%Times are computed cumulatively from the root. If B=0, times are not
%adjusted

GlobalSwitches;

if ~isempty(s(i).child);
   c1=s(i).child(1);
   c2=s(i).child(2);
   s(c1).parent=i;
   if B, s(c1).time=s(c1).time+s(i).time; end
   s(c2).parent=i;
   if B, s(c2).time=s(c2).time+s(i).time; end
   s=FinishBuilding(s,c1,B);
   s=FinishBuilding(s,c2,B);
   s(i).type=ANST;
else
   s(i).type=LEAF;
end


 
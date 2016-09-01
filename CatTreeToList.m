function [s ListOfCats]=CatTreeToList(scat,s)
% Takes a catastrophe tree (as given by the nexus files) scat and a
% corresponding catastrophe-free tree s, and adds the catastrophes to s.

global ANST LEAF

N=length(s);

ListOfCats=zeros(1,N);

for i=1:N
    s(i).cat=[];
    if s(i).type==ANST
        ListOfCats(i)=round(scat(scat(i).parent).time-scat(i).time);
    elseif s(i).type==LEAF
        ListOfCats(i)=round(scat(scat(i).parent).time-scat(i).time-.1);
    end
end

for i=find(ListOfCats)
    s(i).cat=[1/(ListOfCats(i)+1):1/(ListOfCats(i)+1):1-1/(ListOfCats(i)+1)].*(s(s(i).parent).time-s(i).time)+s(i).time;
end
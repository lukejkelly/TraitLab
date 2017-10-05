function [s, clade]=consensuscat(filename,p,OUTGROUP,step,first)
% consensuscat(filename,p,OUTGROUP,step,first)
%draw a consensus trees with catastrophes

global ANST LEAF MCMCCAT
MCMCCAT=1;

if nargin==5
    [s,clade]=consensus([filename,'.nex'],p,OUTGROUP,step,first);
    scat=consensus([filename,'cat.nex'],p,OUTGROUP,step,first);
elseif nargin==4
    [s,clade]=consensus([filename,'.nex'],p,OUTGROUP,step);
    scat=consensus([filename,'cat.nex'],p,OUTGROUP,step);
elseif nargin==3
    [s,clade]=consensus([filename,'.nex'],p,OUTGROUP);
    scat=consensus([filename,'cat.nex'],p,OUTGROUP);
else
    [s,clade]=consensus([filename,'.nex'],p);
    scat=consensus([filename,'cat.nex'],p);
end
N=length(s);

ListOfCats=zeros(1,N);

for i=1:N
    s(i).cat=[];
    if s(i).type==ANST
        ListOfCats(i)=scat(scat(i).parent).time-scat(i).time;
    elseif s(i).type==LEAF
        ListOfCats(i)=scat(scat(i).parent).time-scat(i).time-.1;
    end
end

ListOfCats=round(ListOfCats);

for i=find(ListOfCats)
    s(i).cat=[1/(ListOfCats(i)+1):1/(ListOfCats(i)+1):1-1/(ListOfCats(i)+1)].*(s(s(i).parent).time-s(i).time)+s(i).time;
end


draw(s,1,3,'Consensus Tree');
        


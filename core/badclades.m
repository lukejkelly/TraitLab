function [bad,verybad,verybadn]=badclades(clade,lang)


nc=size(clade,2);
verybad=[];
verybadn=[];
for k=1:nc
    clang=clade{k}.language;
    ncl=size(clang,2);
    tot=ncl;
    bad(k).langs=[];
    bad(k).index=0;
    for j=1:ncl
        if ~any(strcmp(clang{j},lang))
            tot=tot-1;
            bad(k).langs=[bad(k).langs,{clang{j}}];
            bad(k).index=k;
        end
    end
    if tot==0
        verybad=[verybad,{clade{k}.name}];
        verybadn=[verybadn,k];
    end
end
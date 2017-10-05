function ShowCognates(CF,L,s,cognate)

GlobalSwitches;

for k=1:L, draw(s,CF,COGS,'Cognate subtrees',k,cognate{k}); pause(1); end


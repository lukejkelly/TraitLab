function s = treeclades(s,c)

% function s = treeclades(s,c)
%
% write the clade data (c) into the tree (s) leaves
% GKN 28/10

global LEAF

n=size(c,2);
leaves=find([s.type]==LEAF);
for lf=leaves, s(lf).clade=[]; end
languages={s(leaves).Name}; 

for k=1:n
    [lang,InClade,InLangs]=intersect(c{k}.language,languages);
    nodes=leaves(InLangs);
    for i=nodes
        s(i).clade=[s(i).clade,k];
        if size(InLangs,1)==1 & ~isempty(c{k}.rootrange) % Luke 28/11/2014 this changed from size(InLangs,2)==1 & ~isempty(c{k}.rootrange) 
            if size(c{k}.rootrange,2)==2
                s(i).leaf_has_timerange=1;
                s(i).timerange=c{k}.rootrange;
            else
                disp('error in treeclades - leaf time range should have two elements > and <');
            end
        end
    end
end

outleaves=find(~[s(leaves).leaf_has_timerange]);
for k=leaves(outleaves)
    s(k).timerange=[s(k).time,s(k).time];
end

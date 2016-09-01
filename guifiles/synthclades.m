function  clade = synthclades(s,numclade,originrootboth,accuracy);
% clade = synthclades(tree,numclade,originrootboth,accuracy);
global LEAF ANST

clade = {pop('clade')};
clade = repmat(clade,1,numclade);

% get all ancestral nodes
eligible = find([s.type]==ANST);

for i = 1:numclade
    clade{i}.name = ['Clade_' num2str(i)];
    % choose node
    chosen = ceil(rand * length(eligible));
    % find leaves below that node
    leaves = below(s,eligible(chosen));
    leaves = leaves(find([s(leaves).type] == LEAF));
    % get names off leaves and put then in the clade definition
    clade{i}.language = {s(leaves).Name};
    switch originrootboth
    case 1
        doroot = 0;
    case 2
        doroot = 1;
    case 3
        doroot = rand < 0.5;
    end
    if doroot
        clade{i}.rootrange = [s(eligible(chosen)).time s(eligible(chosen)).time] .* [1-accuracy 1+accuracy];
    else
        p = s(eligible(chosen)).parent;
        clade{i}.adamrange = [s(p).time s(p).time].* [1-accuracy 1+accuracy] ;
    end
    % update list of eligible nodes
    eligible(chosen) = [];
end
    
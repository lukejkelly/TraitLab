function ok = checkcladetreematch(s,clade)
global LEAF
ok = 1;
% names of languages from leaves
l = find([s.type] == LEAF);
langname = {s(l).Name};
notfound = {};
for i = 1:length(clade)
    for j = 1:length(clade{i}.language)
        if ~any(strcmp(clade{i}.language(j),langname))
            notfound = [notfound,{clade{i}.name},clade{i}.language(j)];
            ok = 0;
        end
    end
end
        
if ~ok
    disp('Problem with building initial tree: language names in clade block don''t match those in data block')
    disp('Check the following language names in the given clades:')
    disp('(note that all names are case sensitive)')
    disp(sprintf('clade %s: %s\n',notfound{:}))
end
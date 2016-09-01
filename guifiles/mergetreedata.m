function [s1,ok] = mergetreedata(s,data,language);

% [s1,ok] = mergetreedata(s,data,language) checks that the data matrix and the 
% tree s1 are compatible and if so, inserts the data into the appropriate leaves.
% The leaves in the tree must be named with the names in the cell array language
% or in the case that language is empty the languages are assumed to be 
% named 1:N down the rows of the data matrix
% ok == 0 if there was an error and the unaltered tree s is returned

% modified 19/1/05 weakened failure conditions to allow
% data to be a superset of leaf names

global LEAF 

ok = 1;
s1=s;
% check we have non-empty tree
if isempty(s)
    disp('Error in mergetreedata: empty tree passed in')
    ok = 0;
end

if isempty(data)
    disp('Warning from mergetreedata: empty data written onto tree');
    disp('This would be appropriate if you wish to sample the prior');
end

% check that data sizes match
if ok
    leaves = find([s.type] == LEAF);
    ns = length(leaves);  
    if (ns > size(data,1)) | (ns < size(data,1) & isempty(language))
        disp('Error in mergetreedata: taxon <-> leaf mapping not defined')
      %  keyboard;
        ok = 0;
        return
    end
    % if ns = size(data,1) but language names are not given, make cell
    % array of languages for names and assume rows map to leaves by
    % numbering
    if isempty(language)
        language = sprintf('%4.0f',1:ns);
        language = cellstr(reshape(language,ns,4));
    end        
    language = strrep(language,' ' ,'');
end

% compare names of leaf nodes and languages 
if ok
    % get leaf node names
    leafname = {s(leaves).Name};
    % remove any spaces
    leafname = strrep(leafname,' ','');
    %sort leafname and languages
    [sortleaf,sortleafind] = sort(leafname);
    [sortlang,sortlangind] = sort(language);
    %sortleaf = reshape(sortleaf,size(sortlang));
    sortleaf = reshape(sortleaf,[ns,1]);
end

% everything seems ok - stick data in the right place
% if ok
%     datacell=num2cell(data(sortlangind(sortleafind),:),2);
%    [s1(leaves).dat]=deal(datacell{:});
% end

if ok
    for j = 1:ns
        dat_row=strmatch(leafname(j),language,'exact');
        if isempty(dat_row)
            disp(['Error in mergetreedata: no taxon data found for tree leaf ',leafname(j)])
            keyboard;
            ok = 0;
        end
        s1(leaves(j)).dat = data(dat_row,:);
    end
end

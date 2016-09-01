function s=dotreeblock(fid);

% need to process commands TRANSLATE and TREE
translated=0;
token=upper(gettoken(fid));
while ~feof(fid) & ~strcmp(token,'END') & ~strcmp(token,'ENDBLOCK')
    switch token
        
    case 'TRANSLATE'
        symbolused={};
        symbolmeaning={};
        token = gettoken(fid);
        while ~feof(fid) & token~=';'
            symbolused{end+1}=token;
            symbolmeaning{end+1}=gettoken(fid);
            token = gettoken(fid);
            if token == ','
                token = gettoken(fid);
            end
        end
            translated=1;
    case 'TREE'
        % need to skip over tree name and possible * symbol
        token = gettoken(fid);
        while ~feof(fid) & token ~= '='
            token=gettoken(fid);
        end
        % at beginning of tree - read it in
        s=readnexustree(fid);
    otherwise
        skipcommand(fid);
    end       
    token=upper(gettoken(fid));
end

if translated
    % need to replace symbolused with symbolmeaning
    for i=1:length(symbolused)
       q=find(strcmpi(symbolused(i),{s.Name})==1);
       if q>0
           s(q).Name=char(symbolmeaning(i));
       else
           disp('Problem in dotreeblock - translation may not have worked.')
       end
   end
end


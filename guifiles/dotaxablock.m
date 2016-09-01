function [langnames,numlangs] = dotaxablock(fid);

numlangs = 0;
langnames = {};
numread = 0;
% want Dimensions command and Taxlabels
token=upper(gettoken(fid));
while ~feof(fid) & ~strcmp(token,'END') & ~strcmp(token,'ENDBLOCK')
    switch token
    case 'DIMENSIONS'
        while ~feof(fid) & token~=';' 
            token=upper(gettoken(fid));
            switch token
                % if the the token is a recognised command pay attention
                % otherwise just keep reading
            case 'NTAX'
                if gettoken(fid)~='='
                    disp('Error in reading TAXA block. NTAX could not be read');
                else
                    numlangs=str2double(gettoken(fid,'+-'));
                end                    
            end
        end
        token=upper(gettoken(fid));
    case 'TAXLABELS'
        if numlangs == 0 
            disp('Error in reading TAXA block.  TAXLABELS command reached before DIMENSIONS defined')
        else
            token=gettoken(fid);
            while ~feof(fid) & token~=';'
                numread = numread+1;
                if numread > numlangs
                    disp(sprintf('Error in reading TAXA block.  Found more than NTAX = %1.0f taxa labels',numlangs));
                else
                    langnames{numread} = token;
                end
                token = gettoken(fid);
            end
            if numlangs ~= length(langnames)
                disp(sprintf('Error in reading TAXA block.  Found only %1.0f taxlabels, expected NTAX = %1.0f',numread,numlangs));
            end
        end
        token=upper(gettoken(fid));
    otherwise
        skipcommand(fid);
        token=upper(gettoken(fid));
    end
end


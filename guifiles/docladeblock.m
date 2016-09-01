function clade = docladeblock(fid);

clade = {};
n = 0;
badclade =[];
% want CLADE command
token=upper(gettoken(fid));
while ~feof(fid) & ~strcmp(token,'END') & ~strcmp(token,'ENDBLOCK')
    switch token
    case 'CLADE'
        n = n+1;
        clade{n} = pop('clade');
        while ~feof(fid) & token~=';' 
            token=upper(gettoken(fid));
            switch token
                % if the the token is a recognised command pay attention
                % otherwise just keep reading
            case 'NAME'
                if gettoken(fid)~='='
                    disp(sprintf('Error in reading CLADES block. NAME of clade %1.0f could not be read',n));
                    badclade = [badclade n];
                else
                    clade{n}.name = gettoken(fid);
                end                    
            case 'ROOTMIN'
                clade{n}.rootrange(1)  = getvalue(fid,sprintf('ROOTMIN of clade %1.0f',n),'CLADES',-1);
            case 'ROOTMAX'
                clade{n}.rootrange(2) = getvalue(fid,sprintf('ROOTMIN of clade %1.0f',n),'CLADES',-1);
            case 'ORIGINATEMIN'
                clade{n}.adamrange(1) = getvalue(fid,sprintf('ROOTMIN of clade %1.0f',n),'CLADES',-1);
            case 'ORIGINATEMAX'
                clade{n}.adamrange(2) = getvalue(fid,sprintf('ROOTMIN of clade %1.0f',n),'CLADES',-1);
            case 'TAXA'
                m = 0;
                if gettoken(fid)~='='
                    disp(sprintf('Error in reading CLADES block. TAXA for clade %1.0f could not be read',n));
                else
                    token = gettoken(fid);
                    while ~feof(fid) & token~=';' 
                        m = m+1;
                        clade{n}.language{m} = token;
                        token = gettoken(fid); 
                    end
                    % get rid of delimiters
                    delims = strmatch(',',clade{n}.language);
                    weirdchars = strmatch('=',clade{n}.language);
                    if ~isempty(weirdchars)
                        disp(sprintf('Error in reading CLADES block. Possible missed semicolon in clade %1.0f.  Clade ignored',n));
                        badclade = [badclade n];
                    end
                    good = 1:length(clade{n}.language);
                    good(delims)= 0;
                    good = find(good);
                    clade{n}.language = clade{n}.language(good);
                end                    
            end
        end 
        token=upper(gettoken(fid));
        % check that clade is ok 
        if ~isempty(clade{n}.rootrange)
            if any(clade{n}.rootrange < 0)
               %error in reading
              badclade = [badclade n];
            end
            if clade{n}.rootrange(2) == 0
                % max not read, set to infinity
                clade{n}.rootrange(2) = inf;
            end                
            if clade{n}.rootrange(1)>clade{n}.rootrange(2)
                disp(sprintf('Clade %1.0f has ROOTMIN > ROOTMAX.  Clade ignored',n));
                badclade = [badclade n];
            end
        end
        if ~isempty(clade{n}.adamrange)
            if any(clade{n}.adamrange < 0)
                %error in reading
                badclade = [badclade n];
            end
            if clade{n}.adamrange(2) == 0
                % max not read, set to infinity
                clade{n}.adamrange(2) = inf;
            end                
            if clade{n}.adamrange(1)>clade{n}.adamrange(2)
                disp(sprintf('Clade %1.0f has ORIGINATEMIN > ORIGINATEMAX.  Clade ignored',n));
                badclade = [badclade n];
            end
        end
        if isempty(clade{n}.language)
            disp(sprintf('Clade %1.0f has no languages in it.  Clade ignored',n));
        end        
    otherwise
        skipcommand(fid);
        token=upper(gettoken(fid));
    end
end

% get rid of bad clades
good = 1:n;                    
good(badclade)= 0;
good = find(good);
clade = clade(good);


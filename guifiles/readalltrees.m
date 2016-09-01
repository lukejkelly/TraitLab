function alltrees = readalltrees(filename);

alltrees={};
ok =1;
    
% try to open the tree file
try
    if ok
        alltrees=textread([filename],'%s','delimiter',';','endofline',';','bufsize',1e6);
        if isempty(alltrees)
            disp(['Error in readoutput: ' filename ' contains no information']) 
            ok = 0;
        end
    end
catch
    disp(lasterr);
    %if ~isempty(strfind(lasterr,'Buffer'))
    %    disp('To avoid a buffer overflow, try copying the trees you want into a ')
    %    disp('nexus file without including extraneous information such as the Data Block')
    %end
    ok = 0;
end

if ok
    % clear any whitespace from trees
    alltrees = strrep(alltrees,sprintf('\n'),'');
    alltrees = strrep(alltrees,sprintf('\r'),''); %Emailing turns \n into \r
    alltrees = strrep(alltrees,sprintf('\t'),'');
    alltrees = strrep(alltrees,sprintf(' '),'');
    % get lines starting with a comment
    allcomments = alltrees(strmatch('[',alltrees));
    if ~isempty(allcomments)
        %attempt to remove comment
        for i=1:length(allcomments)
            finish = min(strfind(allcomments{i},']'));
            if ~isempty(finish) & length(allcomments(i))>finish
                allcomments{i}=allcomments{i}((finish+1):end);
            end
        end
        % get lines that now start with tree command
        allcomments = allcomments(strmatch('tree',lower(allcomments)));
    end
    % get lines with starting with tree command
    alltrees = alltrees(strmatch('tree',lower(alltrees)));
    alltrees = [alltrees;allcomments]; 
    if isempty(alltrees) 
        disp(['No trees found in ' filename]);
        disp('If there are really trees in the file, try removing any comments in the')
        disp('tree block preceding any TREE commands (comments such as [&R] or [&U] may remain)')
        ok = 0;
    else
        ntree = length(alltrees);
    end
end

if ok
    good = ones(ntree,1);
    % trim front bits off tree strings
    for i=1:ntree
        start = min(strfind(alltrees{i},'('));
        finish = max(strfind(alltrees{i},')'));
        hasequals = strfind(alltrees{i},'=');
        if isempty(hasequals)
                      disp(['No ''='' found in tree ' num2str(i) ' in ' filename '.  Tree ignored']);
            good(i) = 0;
        end 
        if isempty(start) & good(i) == 1
            disp(['No ''('' found in tree ' num2str(i) ' in ' filename '.  Tree ignored']);
            good(i) = 0;
        end 
        if isempty(finish)& good(i) == 1
            disp(['No '')'' found in tree ' num2str(i) ' in ' filename '.  Tree ignored' ]);
            good(i)=0;
        end 
        if start > finish& good(i) == 1
            disp(['First ''('' comes after last '')'' in tree ' num2str(i) ' in ' filename '.  Tree ignored']);
            good(i)=0;
        end 
        if good(i)
            alltrees{i}=alltrees{i}(start:finish);
        end
    end
    alltrees = alltrees(find(good));
end


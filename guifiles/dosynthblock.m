function true = dosynthblock(fid);

true = pop('true');
%set default values for reading in data
[true.mu,true.br,true.lambda,true.theta,true.vocabsize] = deal(-1);

% want PARAMETERS command
token=upper(gettoken(fid));
while ~feof(fid) & ~strcmp(token,'END') & ~strcmp(token,'ENDBLOCK')
    switch token
    case 'PARAMETERS'
        while ~feof(fid) & token~=';' 
            token=upper(gettoken(fid));
            switch token
                % if the the token is a recognised command pay attention
                % otherwise just keep reading
            case 'MU'
                true.mu = getvalue(fid,'Mu','SYNTHESIZE',-1);
            case 'BORROWRATE'
                true.br = getvalue(fid,'Borrowrate','SYNTHESIZE',-1);
            case 'LAMBDA'
                true.lambda = getvalue(fid,'Lambda','SYNTHESIZE',-1);
            case 'THETA'
                true.theta = getvalue(fid,'Theta','SYNTHESIZE',-1);
            case 'VOCABSIZE'
                true.vocabsize = getvalue(fid,'Vocabsize','SYNTHESIZE',-1);
            case 'P'
                true.p = getvalue(fid,'p','SYNTHESIZE',-1);
            case 'BETA' % LUKE 04/09/2016
                true.beta = getvalue(fid, 'beta', 'SYNTHESIZE', -1);
            end
        end
        token=upper(gettoken(fid));
    otherwise
        skipcommand(fid);
        token=upper(gettoken(fid));
    end
    % check that all the statistics are numeric
    toread = {'Mu','Borrowrate','Lambda','Vocabsize','Theta','p', 'beta'};
    notAN = find(isnan([true.mu,true.br,true.lambda,true.vocabsize,true.theta, true.beta]));
    if  ~isempty(notAN);
        disp(sprintf('Error in reading SYNTHESIZE block: %s is not numeric\n', toread{notAN} ));
    end
    % check that all the necessary statistics were found
    notdone = find( [true.mu,true.br,true.lambda,true.vocabsize,true.p, true.beta] < 0 );
    if  ~isempty(notdone);
        disp(sprintf('Error in reading SYNTHESIZE block: %s not found\n', toread{notdone} ));
    end
    % if theta was not found, make it 0
    if true.theta<0
        true.theta=0;
    end
end
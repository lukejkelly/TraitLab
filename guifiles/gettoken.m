function token=gettoken(fid,ignorepunc);

% gettoken(fid,ignorepunc);
% Returns the next token in the nexus file identified by fid
% The current position of the file after running the 
% function is one character after the returned token
% Use the ignorepunc input when reading in a token that is allowed to
% contain otherwise restricted symbols.  ignorepunc is a string of std punctuation
% symbols that are to be treated as plain characters
% example:  to read the numbers 2.441e-10 or 3e+5 call gettoken(fid,'+-')

punc='(){}/\,;:=*''"`+-<>'; %punctuation characters
if nargin==2
    for i=1:length(ignorepunc)
        punc=strrep(punc,ignorepunc(i),'');
    end
end
token='';
charsize=1; %number of bytes in a character
currchar=nextchar(fid);

% need to skip preceding whitespace and comments
while ~feof(fid) & (isspace(currchar) | currchar=='[')
    if currchar=='['  % at start of comment so skip to end
        skipcomment(fid);
        currchar=nextchar(fid);
    else % in patch of whitespace - try next character
        currchar=nextchar(fid);
    end
end

% at start of token - need to check whether it is a quoted token or not

if ~isempty(currchar) & currchar==''''
    % quoted token
    while ~feof(fid)
        currchar=nextchar(fid);
        if isempty(currchar) | currchar~=''''
            % plain character in quote - add to string
            token = [token currchar];
        else
            % quote mark - may be end or may be double
            currchar=nextchar(fid);
            if ~isempty(currchar) & currchar==''''
                % double - add and continue
                token = [token currchar];
            else
                % single - we are done
                % need to back up one character
                fseek(fid,-charsize,'cof');              
                break
            end
        end
    end
elseif ~isempty(currchar) & any(currchar==[punc])
    % a punctuation token
    token = currchar;
else
    % a plain vanilla token - read it all
    %stop if we run into whitespace or punctuation
    while ~feof(fid) & ~isspace(currchar)
        if  any(currchar == punc)
            % have run into punctuation - need to back up one step in file
            fseek(fid,-charsize,'cof');
            %exit loop
            break
        else
            % legitimate part of word - add to end
            token = [token currchar];
        end
        currchar = nextchar(fid);
    end
end

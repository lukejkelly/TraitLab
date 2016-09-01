function skipcomment(fid);

% Advances the current position in the nexus file identified by fid
% to the end of the comment.  
% Assumes that the current position of fid is the beginning of a comment
% Note that nested comments have been catered for

currchar=nextchar(fid);
while ~feof(fid) & currchar~=']' 
    currchar=nextchar(fid);
    if ~isempty(currchar) & currchar=='['    
        skipcomment(fid);
    end
end


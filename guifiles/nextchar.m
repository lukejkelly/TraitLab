function c=nextchar(fid)

% returns the character at the current file position of
% the file identified by fid and advances the current position
% to the next character

c=fscanf(fid,'%c',1);
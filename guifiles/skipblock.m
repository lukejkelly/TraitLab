function skipblock(fid);

token=gettoken(fid);        
while ~feof(fid) & ~strcmpi(token,'END') & ~strcmpi(token,'ENDBLOCK')
    token=gettoken(fid);
end
skipcommand(fid);

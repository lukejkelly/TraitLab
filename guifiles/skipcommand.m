function skipcommand(fid);

token=gettoken(fid);
while ~feof(fid) & token ~= ';' 
    token=gettoken(fid);
end

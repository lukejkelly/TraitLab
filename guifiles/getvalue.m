function value = getvalue(fid,valuename,blockname,default)
value = default;
if gettoken(fid)~='='
    disp(sprintf('Error in reading %s block. %s could not be read',blockname,valuename));
else
    value=str2double(gettoken(fid,'+-'));
end                    

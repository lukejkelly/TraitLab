function a = readrunfile(fname)

fid = fopen(fname);
if fid > 0
a = textscan(fid,'%s%[^\n]%c','delimiter','=','commentStyle','%');
fclose(fid);
a{1} = deblank(a{1});
a = a(1:2);
else
    error('Couldn''t find file %s.  Make sure it is in the correct directory',fname)
end
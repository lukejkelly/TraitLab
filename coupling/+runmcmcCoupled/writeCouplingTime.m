function writeCouplingTime(handles, t)
    fn = [handles.output.path, handles.output.file, '.tau'];
    fid = fopen(fn, 'w+');
    fprintf(fid, '%d\n', t);
    fclose(fid);
end

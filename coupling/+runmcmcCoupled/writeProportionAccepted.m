function writeProportionAccepted(pa, outputPath, outputFile, t)
    fn = [outputPath, outputFile, '.pa'];
    if nargin == 3
        fid = fopen(fn, 'w');
        fprintf(fid, '%s', 't');
        fprintf(fid, ',%d', 1:length(pa));
    else
        fid = fopen(fn, 'a');
        fprintf(fid, '%d', t);
        fprintf(fid, ',%.04f', pa);
    end
    fprintf(fid, '\n');
    fclose(fid);
end

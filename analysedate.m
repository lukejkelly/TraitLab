function res=analysedate(filename,taxon)

[out,ok]=readoutput(filename);
if ~ok
    disp('Error: cannot read output');
    keyboard;
end

Nsamp=length(out.trees);

for i=1:Nsamp
    s=rnextree(out.trees{i});
    found=0;
    cur=0;
    while ~found
        cur=cur+1;
        if strcmp(s(cur).Name,taxon), found=1; end
    end
     
    res(i)=s(cur).time;
end
function res=analysedate2(filename,taxon1,taxon2)

[out,ok]=readoutput(filename);
if ~ok
    disp('Error: cannot read output');
    keyboard;
end

Nsamp=length(out.trees);

for i=1:Nsamp
    s=rnextree(out.trees{i});
    found=0;
    t1=0;
    while ~found
        t1=t1+1;
        if strcmp(s(t1).Name,taxon1), found=1; end
    end
    
    found2=0;
    t2=0;
    while ~found2
        t2=t2+1;
        if strcmp(s(t2).Name,taxon2), found2=1; end
    end
     
    p=mrca([t1,t2],s);
    res(i)=s(p).time;
    
    
%     if ~strcmp(s(s(s(cur).parent).child(3-s(cur).sibling)).Name,taxon2)
%         res(i)=NaN;
%     else
%         res(i)=s(s(cur).parent).time;
%     end
end
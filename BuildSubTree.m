function s=BuildSubTree(s,mask,language)

global LEAF OTHER ROOT ADAM

[s.mark]=deal(0);

%find leaves to be deleted
LEAVES=find([s.type]==LEAF);
NS=length(LEAVES);
for j=LEAVES
    leafname = strrep(s(j).Name,' ','');
    for k=mask
        masklang = strrep(language{k},' ','');
        if isequal(leafname,masklang)
            s(j).mark=1;
            break;
        end
    end
end  

d=find([s.mark]==1);
if (length(d)+2) >= size(s,2)/2
    disp('outtree from BuildSubTree() would have less than two leaves');
    disp('BuildSubTree() returning input tree unchanged');
    return
end
while ~isempty(d)
    
    so=s;
    
    j=d(1);
    s=swap(s,j,2*NS);
    j=2*NS;
    sj=s(j);
    
    %result=checktree(s,NS);
    %if ~isempty(result), disp('error j'); keyboard;pause; end
    
    p=sj.parent;
    s=swap(s,p,2*NS-1);
    p=2*NS-1;
    sp=s(p);
    
    %result=checktree(s,NS);
    %if ~isempty(result), disp('error p'); keyboard;pause; end
    
    soo=s;
    
    pp=s(p).parent;
    poc=s(p).child(OTHER(s(j).sibling));
    s(pp).child(s(p).sibling)=poc;
    s(poc).parent=pp;
    s(poc).sibling=s(p).sibling;    
    if s(p).type==ROOT, 
        if s(poc).type==LEAF
            disp('problem');
        else
            s(poc).type=ROOT;
        end
    end    
    s((2*NS-1):(2*NS))=[];
    NS=NS-1;
    
    %result=checktree(s,NS);
    %if ~isempty(result), disp('error pp'); keyboard;pause; end
        
    a=find([s.type]==ADAM);
    s=swap(s,a,2*NS);
    
    %result=checktree(s,NS);
    %if ~isempty(result), disp('error a'); keyboard;pause; end
    
    d=find([s.mark]==1);
end
result=checktree(s,NS);
if ~isempty(result), disp('BuildSubTree() built a bad tree'); keyboard; end
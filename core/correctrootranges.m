function clade=correctrootranges(clade,i)

%correctrootranges make clade{:}.rootrange calibration intervals consistent
%
%% function clade=correctrootranges(clade,i)
%
%makes the clade constraints consistent and computes the 
%clade constraint for each clade GKN 31/01/2008
%
%Input
%clade, a cell array of TraitLab clades clade{1:end}
%
%Output
%clade, a cell array of TraitLab clades clade{1:end}
%with clades for each leaf and the clade of all taxa.
%Every clade now has a rootrange, and rootranges are clipped to
%be as short as possible in the given heirarchy of clades
%
%% Example
%[ajunk,bjunk,cjunk,clade] = nexus2stype('../a400/a400.nex');
%DATA found
%Reading MATRIX
%CLADES found
% found .... and ignored
%s = Clade2Tree(clade);
%draw(s,2,0,'',[],[],clade); 
%
%See also consensus, draw, nexus2stype

%% recursive - starts top down and ends bottom up
for m=clade{i}.children
    
%% Bring down information from above
    if (isfield(clade{i},'rootrange') && length(clade{i}.rootrange)==2)
        if  (isfield(clade{m},'rootrange') && length(clade{m}.rootrange)==2)
            if clade{m}.rootrange(2)>clade{i}.rootrange(2), clade{m}.rootrange(2)=clade{i}.rootrange(2); end
        else
            clade{m}.rootrange=[0,clade{i}.rootrange(2)];
        end
    end
    
%% Ensure the tree below is up to date
    if length(clade{m}.language)>1
        clade=correctrootranges(clade,m);
    else
        if (~isfield(clade{m},'rootrange') || length(clade{m}.rootrange)~=2), clade{m}.rootrange=[0,inf]; end
    end

%% Bring up information from the tree below
    if (isfield(clade{i},'rootrange') && length(clade{i}.rootrange)==2)
        if clade{m}.rootrange(1)>clade{i}.rootrange(1), clade{i}.rootrange(1)=clade{m}.rootrange(1); end
    else
        clade{i}.rootrange=[clade{m}.rootrange(1),inf];
    end
    
end
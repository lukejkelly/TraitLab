function s = Clade2Tree(clade)
%Clade2Tree convert a clade list into a (consensus) tree
%
%function s = Clade2Tree(clade)
%
%Computes a (multifurcating) tree representing the hierarchy of
%clades/splits. edge length (above clade{i}) is given in clade{i}.length,
%clade{i}.support gives the probability for clade{i}. GKN 10/07/2007
%
%Input
%clade, a cell array of TraitLab clades clade{1:end}
%
%Output
%s, a TraitLab tree structure.  clade{i} generates node s(i).
%s is a binary tree, with edges of length "eps" to imitate
%multifurcation for drawing purposes (this would be OK under
%non-cladagenic simulation). For ancestral nodes the s.Name
%field is used to hold the percent support. Values under 95%
%only are shown.
%
%%Example 1: used by consensus to make consensus trees
%
%%Example 2: Visualising arbitrary clade structures. Suppose
%%IEcnsrvSwd100.nex is a nexus file in the current path containing a
%%TraitLab CLADE block
%
% [ajunk,bjunk,cjunk,clade] = nexus2stype('IEcnsrvSwd100.nex');
% s = Clade2Tree(clade);
% draw(s,2,3,'Trimmed Consensus Tree'); %figure(2), CONC=3 draw as CT
%
%See also consensus, draw, nexus2stype

global ADAM ROOT


%% get all taxa
nc=length(clade);
langlot=[];
NOLENGTH=0;
UPBOUNDCLADE=0;
for k=1:nc
    langlot=[langlot, clade{k}.language];
    if ~isfield(clade{k},'length'), NOLENGTH=1; clade{k}.length=1; end
    if ~isfield(clade{k},'support'), clade{k}.support=1; end
    %if ~isfield(clade{k},'children'), clade{k}.children=[]; end;
    %next ignores possible UPBOUNDCLADE constraints from 'adamrange'-ed nodes
    if (isfield(clade{k},'rootrange') && length(clade{k}.rootrange)==2 && ~isinf(clade{k}.rootrange(2))), UPBOUNDCLADE=1; end
end
DO_ROOTRANGE_AGES=(NOLENGTH & UPBOUNDCLADE);

langlot=unique(langlot);
langlot_r=langlot;
nl=length(langlot);

%% find missing clades
topclade=[];
for k=1:nc
    if ischar(clade{k}.language), clade{k}.language={clade{k}.language}; end
    if length(clade{k}.language)==1
        langlot=setdiff(langlot,clade{k}.language);
    end
    if length(clade{k}.language)==nl
        topclade=k;
    end
end

%% fix topclade
if isempty(topclade)
    clade{nc+1}.language=langlot_r;
    clade{nc+1}.name='topclade';
    clade{nc+1}.rootrange=[];
    clade{nc+1}.adamrange=[];
    clade{nc+1}.length=1;
    clade{nc+1}.support=1;
    nc=length(clade);
    topclade=nc;
end

%% fix leaf clades
for i=1:length(langlot)
    clade{nc+i}.language={langlot{i}}; %a cellstr of length one, not a char
    clade{nc+i}.name=langlot{i};
    clade{nc+i}.rootrange=[];
    clade{nc+i}.adamrange=[];
    clade{nc+i}.length=1;
    clade{nc+i}.support=1;
end
nc=length(clade);

%% make clade children-heirarchy
%make the clades into a tree. clade a is a child clade b if the set of taxa of
%clade a are a subset of the set of taxa of clade b. topclade is the root clade
leafclade=[];
for k=1:nc
    clade{k}.children=[];
    if length(clade{k}.language)==1
        leafclade=[leafclade,k];
    else
        for m=1:nc
            if m~=k && isempty(setdiff(clade{m}.language,clade{k}.language))
                ADDTHISONE=1;
                ckc=clade{k}.children;
                for n=1:length(ckc)
                    if isempty(setdiff(clade{ckc(n)}.language,clade{m}.language)), clade{k}.children(clade{k}.children==ckc(n))=[]; end
                    if isempty(setdiff(clade{m}.language,clade{ckc(n)}.language)), ADDTHISONE=0; break; end
                end
                if ADDTHISONE, clade{k}.children(end+1)=m; end
            end
        end
    end
end

%%use rootranges to set times
if DO_ROOTRANGE_AGES
    clade=correctrootranges(clade,topclade);
    if isinf(clade{topclade}.rootrange(2))
        mxrr=0;
        for k=1:nc,
            if ~isinf(clade{k}.rootrange(2)), mxrr=max(mxrr,clade{k}.rootrange(2)); end
        end
    end
end

%% clade2tree
%turn the tree structure in clade{}.children into an s-structure
%write in only the child information (not the parent information)
%write the edge length into the time field, these will be accumulated
%from the root (and then reversed) to make time depth info
node=nc;
for n=1:nc
    s(n)=pop('node');
    if DO_ROOTRANGE_AGES
        if isinf(clade{n}.rootrange(2))
            s(n).time=1.5*mxrr;
        else
            s(n).time=mean(clade{n}.rootrange);
        end
        ctime=s(n).time;
    else
        s(n).time=clade{n}.length;
    end
    if length(clade{n}.language)>1 %at least 2 taxa in this clade, so it is ancestral
        if clade{n}.support<0.995, s(n).Name=num2str(round(100*clade{n}.support)); else, s(n).Name=''; end
        cn=clade{n}.children(1); %1st clade child one goes in as an s-child
        s(n).child=cn;
        i=n; %node i is the 'current' node to which we are adding children.
        for m=2:length(clade{n}.children)
            if m==length(clade{n}.children)
                %last clade-child so drop this in as an s-child of i
                s(i).child(2)=clade{n}.children(m);
            else
                %more than 1 clade children to come, so pop a new s-node
                %called node, hang it in as the second child of i, and give node
                %the next clade child on the list as node's first child
                node=node+1;
                s(node)=pop('node');
                s(node).child=clade{n}.children(m);
                if UPBOUNDCLADE && NOLENGTH, ctime=ctime-ctime*eps; s(node).time=ctime; else, s(node).time=eps; end
                s(node).Name='';
                s(i).child(2)=node;
                i=node;
            end
        end
    else
        %clade n has one taxon so it generates a leaf, its parent was handled above
        s(n).Name=clade{n}.language{1};
    end
end

%% finish as a TraitLab tree
%make is a proper s-tree with an adam-node
s(end+1)=pop('node');
if ~DO_ROOTRANGE_AGES, s(end).time=-realmax; end
s(end).type=ADAM;
s(end).child=topclade;
s(topclade).parent=length(s);

%put in all the rest of the tree info
if ~DO_ROOTRANGE_AGES, s(topclade).time=0; end
s=FinishBuilding(s,topclade,~DO_ROOTRANGE_AGES);
s(topclade).type=ROOT;

%reverse direction of time from the deepest leaf
if ~DO_ROOTRANGE_AGES
    mt=max([s.time]);
    for k=1:length(s), s(k).time=mt-s(k).time; end
end

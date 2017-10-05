function [s,clade]=consensus(filename,p,OUTGROUP,step,first)

% consensus computes the consensus tree of trees in nexus file
%
% function [s,clade]=consensus(filename,p,step,first)
%
% Computes a consensus tree for the trees in the nexus file
% "filename" with threshold "p" (so p=0.5 is majority rule CT).
% Writes to filename.con. For large files (many trees with
% many leaves) you may encounter memory problems. You may wish
% to thin the collection of trees. If "step" is defined, then
% every step-th tree is taken from filename. GKN 10/07/2007
%
% If "first" is defined, all trees between 1 and "first"-1 will be ignored.
% RJR 12/09/2009
%
% Input
% filename, a string giving a path to a nexus file of trees
% p, (optional default=0.5) a real scalar 0.5 <= p <= 1
% OUTGROUP, (optional) a string leafname of the outgroup
% step, (optional) an integer subsample interval, step >= 1
%
% Output
% s, a TraitLab tree structure
% clade, a cell array of clades , clade{1:end}
% (and s is written to filename.con)
%
% %Example, suppose tlout.nex is a file of nexus trees
% [s,clade]=consensus('tlout.nex',0.5); %majority rule CT
% draw(s,1,3,'Consensus Tree'); %figure(1), CONC=3 draw as CT
% [s,clade]=consensus('tlout.nex',0.5,'',2); %take every 2nd tree
% draw(s,2,3,'Consensus Tree from thinned MCMC output');
%
% See also draw, Clade2Tree, consensuscat

%% Globals
%We need the ROOT ANST LEAF ADAM definitions
GlobalSwitches;
GlobalValues;

%% load trees
if ~exist('p','var') || isempty(p),
   p=0.5;
   disp(sprintf('Computing majority rule consensus tree, p=0.5'));
else
   if (p<0.5 || p>1)
       error('p should be in the interval [0.5,1]');
   end
end
disp(sprintf('1 of 5: loading trees'));
trees = readalltrees(filename);
disp(sprintf('Found %d trees in %s',length(trees),filename));

if exist('first','var') && ~isempty(first)
   if (first<1 || first>length(trees)), error('Bad first tree value argument. %d trees in file %s.',length(trees),filename); end
   trees={trees{first:end}};
   if first>1, disp(sprintf('Ignoring first %d trees leaves %d trees',first-1,length(trees))); end
end

if exist('step','var') && ~isempty(step)
   if (step<1 || step>length(trees)), error('bad step value argument'); end
   trees={trees{1:step:end}};
   if step>1, disp(sprintf('Subsampling at step %d leaves %d trees',step,length(trees))); end
end
nt=length(trees);
if nt<2
   error(['Less than 2 trees found in ',filename]);
end

%%allocate memory for the (nt*nn) x (nl+1) binary splits-matrix
t=rnextree(trees{1});
leaves=find([t(:).type]==LEAF);
nl=length(leaves);
nn=(2*nl-1);
leafnames={t(leaves).Name};
if (exist('OUTGROUP','var') && ~isempty(OUTGROUP))
  if ~ischar(OUTGROUP),
      error('OUTGROUP argument is present and not string or []');
  end
  OGi=find(strcmp(OUTGROUP,leafnames));
  if length(OGi)~=1,
       error(sprintf('found %d leaves matching OUTGROUP %s',length(OGi),OUTGROUP));
  end
  useOUTGROUP=1;
else
  useOUTGROUP=0;
end
r=zeros(nt*nn,nl+1);

%% compute all splits of all trees
%nt trees with nl leaves and nn nodes. go through them. each edge in each
%tree has a 1 x nl binary representation marking the split. r is a matrix
%of all these splits. the first column of r give the length of the edge so
%r is (nt*nn) x (nl+1)
disp(sprintf('2 of 5: extract splits from trees \n'));
pc=ceil(nt/100);
for i=1:nt
   if ~mod(i,10*pc), disp(sprintf('processed %4d percent of trees',round(i/pc))); end
   t=rnextree(trees{i});
   root=find([t(:).type]==ROOT);
   %edges() returns a (nn x nl) matrix of splits plus an extra first
   %column of edge lengths
   rdbl=edges(t,root,leafnames);
   if useOUTGROUP
       %standardize the splits for leaf 1 as an outgroup
       for ri=1:nn
           if rdbl(ri,1+OGi)==1,
               ste=sum(rdbl(ri,2:end));
               if ~(ste==1 || ste==nl),
                   rdbl(ri,2:end)=1-rdbl(ri,2:end);
               end;
           end
       end
   end
   r( (((i-1)*nn)+1):(i*nn),:)=rdbl;
end

%% extract list of distinct splits
disp(sprintf('\n3 of 5: form consensus array\n'));
%the memory hungry bit - compute the unique splits
[v,i,j]=unique(r(:,2:end),'rows');
%j(k) says where in v the kth row of r can be found so a(m) is the
%number of times the m'th row of v appeared
[a,b]=hist(j,1:length(i));
els=zeros(1,length(i));
sup=els;
for ii=1:length(i), vels(ii)=mean(r(j==ii,1)); vsup(ii)=sum(j==ii)/nt; end
%v, vels and vsup are the split matrix (rows are splits), mean split length
%and probability of split, w, els and sup are the same thing, thinned so that
%only splits with probability at least p are included
w=v(((a./nt)>p),:);
els=vels((a./nt)>p);
sup=vsup((a./nt)>p);

%% make clades from binary split matrix
%from here on we have the job of converting the split matrix into a tree
%convert the split matrix to clades - notice the extra fields "length"
%and "support".
disp(sprintf('4 of 5: build clades\n'));
clade=cell(1,size(w,1));
for k=1:size(w,1)
   clade{k}.language=leafnames(find(w(k,:)));
   clade{k}.name=['scribble', num2str(k)];
   clade{k}.rootrange=[];
   clade{k}.adamrange=[];
   clade{k}.length=els(k);
   clade{k}.support=sup(k);
end


%% make consensus tree from clade list
disp(sprintf('5 of 5: convert clades to tree\n'));
s = Clade2Tree(clade);


%% save CT
[str,ok] = stype2nexus(s,['Consensus Tree for ',filename],'TREE');
if ok
   fid=fopen(strcat(filename,'.con'),'w');
   if fid>0
       fprintf(fid,str);
       fclose(fid);
   else
       warning(['Could not write consensus tree to ',filename,'.con']);
   end
else
   warning(['Could not write consensus tree to ',filename,'.con']);
end



%% binary splits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=edges(s,i,leafnames)

%edges compute splits of tree s, return as splits x leaves binary matrix

els=s(s(i).parent).time-s(i).time;
if isempty(s(i).child)
   r=zeros(1,length(leafnames));
   i=find(strcmp(s(i).Name,leafnames));
   r(i)=1;
   r=[els,r];
else
   c1=s(i).child(1);
   c2=s(i).child(2);
   r1=edges(s,c1,leafnames);
   r2=edges(s,c2,leafnames);
   r=[[els,(sum([r1(:,2:end);r2(:,2:end)])>eps)];r1;r2];
end



function alltrees = readalltrees(filename);

alltrees={};
ok =1;

% try to open the tree file
try
   if ok
       alltrees=textread([filename],'%s','delimiter',';','endofline',';','bufsize',1e6);
       if isempty(alltrees)
           disp(['Error in readoutput: ' filename ' contains no information'])
           ok = 0;
       end
   end
catch
   disp(lasterr);
   %if ~isempty(strfind(lasterr,'Buffer'))
   %    disp('To avoid a buffer overflow, try copying the trees you want into a ')
   %    disp('nexus file without including extraneous information such as the Data Block')
   %end
   ok = 0;
end

if ok
   % clear any whitespace from trees
   alltrees = strrep(alltrees,sprintf('\n'),'');
   alltrees = strrep(alltrees,sprintf('\r'),''); %Emailing turns \n into \r
   alltrees = strrep(alltrees,sprintf('\t'),'');
   alltrees = strrep(alltrees,sprintf(' '),'');
   % get lines starting with a comment
   allcomments = alltrees(strmatch('[',alltrees));
   if ~isempty(allcomments)
       %attempt to remove comment
       for i=1:length(allcomments)
           finish = min(strfind(allcomments{i},']'));
           if ~isempty(finish) & length(allcomments(i))>finish
               allcomments{i}=allcomments{i}((finish+1):end);
           end
       end
       % get lines that now start with tree command
       allcomments = allcomments(strmatch('tree',lower(allcomments)));
   end
   % get lines with starting with tree command
   alltrees = alltrees(strmatch('tree',lower(alltrees)));
   alltrees = [alltrees;allcomments];
   if isempty(alltrees)
       disp(['No trees found in ' filename]);
       disp('If there are really trees in the file, try removing any comments in the')
       disp('tree block preceding any TREE commands (comments such as [&R] or [&U] may remain)')
       ok = 0;
   else
       ntree = length(alltrees);
   end
end

if ok
   good = ones(ntree,1);
   % trim front bits off tree strings
   for i=1:ntree
       start = min(strfind(alltrees{i},'('));
       finish = max(strfind(alltrees{i},')'));
       hasequals = strfind(alltrees{i},'=');
       if isempty(hasequals)
                     disp(['No ''='' found in tree ' num2str(i) ' in ' filename '.  Tree ignored']);
           good(i) = 0;
       end
       if isempty(start) & good(i) == 1
           disp(['No ''('' found in tree ' num2str(i) ' in ' filename '.  Tree ignored']);
           good(i) = 0;
       end
       if isempty(finish)& good(i) == 1
           disp(['No '')'' found in tree ' num2str(i) ' in ' filename '.  Tree ignored' ]);
           good(i)=0;
       end
       if start > finish& good(i) == 1
           disp(['First ''('' comes after last '')'' in tree ' num2str(i) ' in ' filename '.  Tree ignored']);
           good(i)=0;
       end
       if good(i)
           alltrees{i}=alltrees{i}(start:finish);
       end
   end
   alltrees = alltrees(find(good));
end

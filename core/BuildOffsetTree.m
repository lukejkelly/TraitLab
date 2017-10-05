function [s,roottime] = BuildOffsetTree(lang,time,rootmin,rootmax,rho)
% s = BuildOffsetTree(lang,time,rootmin,rootmax)
% builds a random tree
% with leaves at specified times

global LEAF ANST ROOT ADAM

[time,i] = sort(time);
lang = lang(i);

timestep = unique(time);
ntime = length(timestep);
% see how many languages there are at each time step
for i = 1:ntime
    nlang(i) = sum(time == timestep(i));
end
totallang = cumsum(nlang);

rootmin = max(rootmin,time(end));
if rootmin>rootmax
    error('Inconsistent clade times - check subclade times fall below superclade bounds');
end

% build a coalescent tree introducing leaves when necesssary
% choose coalescing rate so it will probably coalesce around the rootmin time
% if it hasn't coalesced by rootmin time then 
if rootmin>0
rate = 1/rootmin;
end
    
n = nlang(1);
t = time(1);
istep = 1;
% first make leaf nodes
for i = 1:n
    s(i) = TreeNode([],[],[],time(i),lang{i},LEAF);
end
availablechild = [1:n];
nextnode = n+1;
if length(timestep) > istep
    timeofnextstep = timestep(istep+1);
else
    timeofnextstep = rootmin;
end
tinc = 0;

% stop when all languages have been introduced and n==1 
while t < rootmin
    % get a time to the next possible coalescence
    if n>1; tinc = -log(rand)/(n*(n-1)*rate); end
    % check that time is before the next bunch of languages 
    % and that there are still 2 or more lineages to coalesce
    if t+tinc > timeofnextstep | n==1
        % gone past next step
        % if next step is not the rootmin
        if length(timestep) > istep
            %introduce more languages
            istep = istep + 1;
            t = timestep(istep);
            n = nlang(istep)+n;
            % get names of languages right
            j = totallang(istep - 1) +1;
            % make nodes for the new languages
            for i = nextnode:nextnode+nlang(istep)-1
                s(i) = TreeNode([],[],[],t,lang{j},LEAF);
                j = j+1;
            end
            availablechild = [availablechild nextnode:nextnode+nlang(istep)-1];
            nextnode = nextnode+nlang(istep);
            % set timeofnextstep
            if length(timestep) > istep
                timeofnextstep = timestep(istep+1);
            else
                timeofnextstep = rootmin;
            end
            
        else
            % hit the rootmin so end it
            t = rootmin;
        end
    else
        % they've coalesced
        t = t+tinc;
        n = n-1;
        % choose two available nodes
        irand = ceil(rand(1,2).*[length(availablechild) length(availablechild)-1]);
        for i = 1:2
            child(i) = availablechild(irand(i));
            availablechild(irand(i)) = [];
        end
        % make new node 
        %%anst nodes dont need names - creating new names is dangerous as
        %%might choose existing leaf name - helpful to have unique names 
        %%unique - was a bug GKN 11/3/11
        %s(nextnode) = TreeNode([],[],child,t,num2str(nextnode),ANST);
        s(nextnode) = TreeNode([],[],child,t,'',ANST);
        % set parent and sibling fields in children
        for i = 1:2
            s(s(nextnode).child(i)).parent = nextnode;
            s(s(nextnode).child(i)).sibling = i;
        end
        % make this new node available as a child node
        availablechild = [availablechild nextnode];
        %[s,availablechild] = makenode(s,availablechild,nextnode,t);
        nextnode = nextnode+1;
    end
end

if n>1
    % force the others to coalesce within the time remaining before rootmax
    % generate uniformly distributed coalescing times in interval [rootmin rootmax]
    ctimes = rand(1,n-1)*(rootmax - rootmin) + rootmin;
    ctimes = sort(ctimes);
    for i = 1:n-1
        [s,availablechild] = makenode(s,availablechild,nextnode,ctimes(i));
        nextnode = nextnode+1;
    end
end

s(end).type = ROOT;
roottime = s(end).time;

% scale if necessary
if roottime < rootmin | roottime > rootmax
    % choose a time uniformly in the allowable root interval
    roottime = rand*(rootmax - rootmin) + rootmin;
    % scale all ANST and ROOT nodes to get root time right
    %check to see whether this tree is just an offset leaf
    if length(lang) > 1
        scalefactor = roottime/s(end).time;
        for i = 1:length(s)
            if any(s(i).type ==[ANST ROOT])
                s(i).time = s(i).time * scalefactor;
            end
        end
    else
        s(1).time = roottime;
    end
end


% make an adam node
s(nextnode) = TreeNode([],[],nextnode-1,realmax,'Adam',ADAM);
s(nextnode-1).parent = nextnode;
s(nextnode-1).sibling = 1;


% %add catastrophies
% cat=zeros(i,1);
% for i=1:totallang(end)
%     if s(i).type==ANST | s(i).type==LEAF
% 	p=s(i).parent;
% 	dt=s(p).time-s(i).time;
% 	cat(i)=poissrnd(rho*dt);
%     end
% end




ok = checktree(s,length(lang));
if ~isempty(ok)
    error('Tree is not legal')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s,availablechild] = makenode(s,availablechild,nextnode,t)

global ANST

% choose two available nodes
irand = ceil(rand(1,2).*[length(availablechild) length(availablechild)-1]);
for i = 1:2
    child(i) = availablechild(irand(i));
    availablechild(irand(i)) = [];
end
% make new node %%GKN 11/3/11 dont think we want a name for a ANST node
% danger of generating names in this way may conflict existing lang name
% was actually happening when working with synth data
% s(nextnode) = TreeNode([],[],child,t,num2str(nextnode),ANST);
s(nextnode) = TreeNode([],[],child,t,'',ANST);
% set parent and sibling fields in children
for i = 1:2
    s(s(nextnode).child(i)).parent = nextnode;
    s(s(nextnode).child(i)).sibling = i;
end
% make this new node available as a child node
availablechild = [availablechild nextnode];

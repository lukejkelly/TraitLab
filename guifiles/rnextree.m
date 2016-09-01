function [s,errmess,ncat]=rnextree(str,zeroleaves,catstr)

% s=rnextree(nexstr) converts the nexus format tree nexstr into
% an s-type tree.  
% It is assumed that str is a WELL FORMED nexus tree
% 
% rnextree(nexstr,1) will set all leaf times to 0
%
% which is assumed to have current position at the first open
% bracket of a tree in nexus format.  Will attempt read tree until
% the first semicolon is reached.
% If nexus tree has any times recorded the times will be left as given 
% so some node times may equal the time of their parent node.
% If no times are recorded, each branch will automatically be 
% assigned length 1 and node times calculated accordingly 

% Modified by RJR 07/11/2008: now also accepts catstr as input. This allows
% to include catastrophes. 
% TODO: check that catstr corresponds to str.

global LEAF ANST ROOT ADAM

if nargin ==1;
    zeroleaves = 0;
end

if nargin < 3
    catstr='';
end


% see how large the tree is
n = 2*sum(str=='(')+2;
% initialise tree of that size
s=TreeNode([],[],[],[],[],[]);
s=repmat(s,1,n);
%convert string to cell array where each entry is a token 
str=strrep(str,'(',' ( ');
str=strrep(str,',',' , ');
str=strrep(str,')',' ) ');
str=strrep(str,':',' : ');
str=strrep(str,';','');
str=strread(str,'%s');

catstr=strrep(catstr,'(',' ( ');
catstr=strrep(catstr,',',' , ');
catstr=strrep(catstr,')',' ) ');
catstr=strrep(catstr,':',' : ');
catstr=strrep(catstr,';','');
catstr=strread(catstr,'%s');

token = 1;

errmess='';
nextnode = 1;
% expected anticipates the next token - tree always starts with (
% note that we use @ to represent a name string
expected = '(';  
nodestack = [];
newsibling=1;
timegiven=0;

%  check first token is bracket
if ~strcmp(str{token}, expected)
    errmess = 'Nexus tree does not start with (';
    disp(errmess);
    pause;
else
    expected='(@';
    % create tree with root as its first node
    currnode = nextnode;    
    nextnode = nextnode+1;
    nodestack = [nodestack currnode];
    
    s(currnode).type = ROOT;
    s(currnode).sibling = newsibling;
    s(currnode).time = 0;
    rootno=currnode;
    % get next token
    token = token+1;
 %   notok=notok+1;tread{notok}=token;
    
end

while token <= length(str);
    if any(str{token}(1)==expected) || (any('@'==expected) && all(str{token}(1)~='():,'))
        % token as expected
        switch str{token}
        case '(' % create new ancestral node            
            % update child and sibling information
            s(currnode).child = [s(currnode).child nextnode];
            newsibling = length(s(currnode).child);
            newparent = currnode;
            currnode = nextnode;
            nextnode = nextnode + 1;
            nodestack = [nodestack currnode]; %#ok<AGROW>
            % create the node
            s(currnode).parent = newparent;
            s(currnode).sibling = newsibling;
            s(currnode).time = 0;
            s(currnode).type = ANST;
                expected='(@';
        case {')',','}  % finished with current node - regress in stack  
            nodestack = nodestack(1:end-1);
            currnode = nodestack(end);
            if str{token}==')'
                expected=',):';
            else % token == ','
                expected='(@';
            end
        case ':' % time data follows
            % Since time data in nexus file is just branch length, simply record
            % this info and make adjustments at end of reading
            token=token+1;
           % gettoken(fid,'+-'); % notok=notok+1;tread{notok}=token;
            
            if token > length(str)
                disp('Unexpected end of nexus string in rnextree - expected time');
                pause
            else
                s(currnode).time = str2num(str{token}); %#ok<ST2NM>
                timegiven = 1;
            end
            expected=',)';
            
            %get number of catastrophes
            if nargin==3
                ncat(currnode)=str2num(catstr{token})-.1; %#ok<ST2NM,AGROW>
            end
            
        otherwise % create leaf node with token as Name
            % update child and sibling information
            s(currnode).child = [s(currnode).child nextnode];
            newsibling = length(s(currnode).child);
            newparent = currnode;
            currnode = nextnode;
            nextnode = nextnode + 1;
            nodestack = [nodestack currnode]; %#ok<AGROW>
            % create the node
            s(currnode).parent = newparent;
            s(currnode).sibling = newsibling;
            s(currnode).time = 0;
            s(currnode).Name = str{token};
            s(currnode).type = LEAF;
            expected=',:)';
            
        end
        token=token+1;%notok=notok+1;tread{notok}=token;
    else
        % unexpected token
        disp(['Error in rnextree - read ' str{token} ' when expected ' expected]);
        keyboard
    end
    
end

% create adam node as last node in tree
adamnode = nextnode;    
s(adamnode).child = rootno;
s(adamnode).Name = 'Adam';
s(adamnode).type = ADAM;
s(adamnode).time = 0;
% name parent in root node
s(rootno).parent=adamnode;

if ~timegiven
    % no times specified - assign arbitrary langth
    [s.time]=deal(1);
end

% adjust times
s=timeconvert(s,rootno);
% find minimum time - ie time furthest from root - and add to all times
timeadj=abs(min([s.time]));
newtime=num2cell([s.time]+timeadj);
[s.time]=deal(newtime{:});
% for i=1:length(s)
%     s(i).time = s(i).time + timeadj;
% end

if ~timegiven && zeroleaves
[s([s.type]==LEAF).time]=deal(0);
end;

s(adamnode).time = realmax;


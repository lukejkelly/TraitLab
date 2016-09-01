function [s,errmess]=readnexustree(fid,zeroleaves)

% s=readnexustree(fid) reads from the file identified by fid
% which is assumed to have current position at the first open
% bracket of a tree in nexus format.  Will attempt read tree until
% the first semicolon is reached.
% If nexus tree has any times recorded the times will be left as given 
% so some node times may equal the time of their parent node.
% If no times are recorded, each branch will automatically be 
% assigned length 1 and node times calculated accordingly 
% readnexustree(fid,1) will set all leaf times to zero
% 
% very similar to rnextree (readnexustree takes a file id, rnextree takes a string)


global LEAF ANST ROOT ADAM

if nargin ==1;
    zeroleaves = 0;
end
errmess='';
nextnode = 1;
% expected anticipates the next token - tree always starts with (
% note that we use @ to represent a name string
expected = '(';  
nodestack = [];
newsibling=1;
timegiven=0;
%notok=0;


%  special case - get first token and check it is bracket
token = gettoken(fid);
%notok=notok+1;tread{notok}=token;
if token ~= expected
    errmess = 'Nexus tree does not start with (';
    disp(errmess);
    pause;
else
    expected='(@';
    % create tree with root as its first node
    currnode = nextnode;    
    nextnode = nextnode+1;
    nodestack = [nodestack currnode];
    s(currnode) = TreeNode([],newsibling,[],0,'',ROOT);
    rootno=currnode;
    % get next token
    token = gettoken(fid);
 %   notok=notok+1;tread{notok}=token;
    
end

while ~feof(fid) && ~strcmp(token,';') %% DW 12 April 2011, token not always a scalar so token~=';' caused error 
    if any(token(1)==expected) || (any('@'==expected) && all(token(1)~='():,'))
        % token as expected
        switch token
        case '(' % create new ancestral node            
            % update child and sibling information
            s(currnode).child = [s(currnode).child nextnode];
            newsibling = length(s(currnode).child);
            newparent = currnode;
            currnode = nextnode;
            nextnode = nextnode + 1;
            nodestack = [nodestack currnode]; %#ok<AGROW>
            % create the node
            s(currnode) = TreeNode(newparent,newsibling,[],0,'',ANST);
                expected='(@';
        case {')',','}  % finished with current node - regress in stack  
            nodestack = nodestack(1:end-1);
            currnode = nodestack(end);
            if token==')'
                expected=',):';
            else % token == ','
                expected='(@';
            end
        case ':' % time data follows
            % Since time data in nexus file is just branch length, simply record
            % this info and make adjustments at end of reading
            token=gettoken(fid,'+-'); % notok=notok+1;tread{notok}=token;
            
            if feof(fid)
                disp('Unexpected end of file in readnexustree - expected time');
                pause
            else
                s(currnode).time = str2num(token); %#ok<ST2NM>
                timegiven = 1;
            end
            expected=',)';
        otherwise % create leaf node with token as Name
            % update child and sibling information
            s(currnode).child = [s(currnode).child nextnode];
            newsibling = length(s(currnode).child);
            newparent = currnode;
            currnode = nextnode;
            nextnode = nextnode + 1;
            nodestack = [nodestack currnode];%#ok<AGROW>
            % create the node
            s(currnode) = TreeNode(newparent,newsibling,[],0,token,LEAF);
            expected=',:)';
            
        end
        token=gettoken(fid);%notok=notok+1;tread{notok}=token;
    else
        % unexpected token
        disp(['Error in readnexustree - read ' token ' when expected ' expected]);
        keyboard
    end
    
end

% create adam node as last node in tree
adamnode = nextnode;    
s(adamnode) = TreeNode([],[],rootno,0,'Adam',ADAM);
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
function [s,content,true,clade]=nexus2stype(filename)

% [s,content,true,clade] = nexus2stype(filename) 
% reads the nexus file filename and extracts relevant information
% to construct an s-type tree structure,s.
% Set the (optional) readdata and readtree inputs to 0  to ignore these blocks
% The default is 1 for both.
%
% If no tree structure is given in the file use the other outputs
% to extract relevant data.  content.array is a numerical matrix with 
% rows = languages and columns = cognate.  
% content.language and content.cognate are cell arrays of strings
% The nth row of content.array corresponds to language content.language(n)
%
% Currently the nexus blocks that it recognises are
% DATA, TAXA, CHARACTERS, TREES, SYNTHESIZE, CLADE 
% All other blocks are ignored
%
% Note that it is necessary to explicitly state values for NTAX and NCHAR
% when nexus file contains a data matrix

global LEAF

if isempty(LEAF)
    GlobalSwitches;
    GlobalValues;
end


ok=1;
s = [];
content = pop('content');
true = pop('true');
clade = {};
treeread=0;

% open file with read permission
[fid,errmess] = fopen(filename,'r');

% check that filename has opened correctly
if fid==-1
    disp(['Unable to open ' filename])
    disp(errmess);
    ok=0;
end

% check that filename is a nexus file
if ok
    token=gettoken(fid);
    ok=strcmpi(token,'#NEXUS');
    if ~ok
        disp([filename ' is not a NEXUS file - missing initial #NEXUS.  Error in nexus2stype']);
    end
end

% read nexus file block by block until end of file
while ok && ~feof(fid)
    % find beginning of block
    token=gettoken(fid);
    while  ~feof(fid) && ~strcmpi(token,'BEGIN') 
        token=gettoken(fid);
    end
    
    % get blockname - needs to be uppercase for switch
    currblock=upper(gettoken(fid));
    skipcommand(fid);
    switch currblock
    case 'DATA'
        disp([currblock ' found']);
        [langnames,content.cognate,content.array] = dodatablock(fid,content.NS);
        if isempty(content.language) 
            content.language = langnames(:);
        end
        content.cognate = content.cognate(:);
        if isempty(content.array)
            ok = 0;
            disp('Matrix not read in DATA block.  Check the matrix is present and the block is correctly formatted.')
        end
    case 'TAXA'
        disp([currblock ' found']);
        [langnames numlangs] = dotaxablock(fid);
        if numlangs > 0 && length(langnames)==numlangs
            content.language = langnames(:);
            content.NS = numlangs;
        end
%     case 'CHARACTERS'
%         disp([currblock ' found']);
%         [langnames,content.cognate,content.array] = dodatablock(fid,content.NS);
%         if isempty(content.language) 
%             content.language = langnames(:);
%         end
%         content.cognate = content.cognate(:);
    case 'TREES'
        disp([currblock ' found']);
        s=dotreeblock(fid);
        treeread=1;
    case 'SYNTHESIZE'
        disp([currblock ' found']);
        true = dosynthblock(fid);
    case 'CLADES'
        disp([currblock ' found']);
        clade = docladeblock(fid);
    case ''
        skipblock(fid);
    otherwise 
        disp([currblock ' found ... and ignored']);
        skipblock(fid);
    end 
end

% close filename
if ok
 [content.NS content.L] =size(content.array);
   if fclose(fid)==-1
        disp(['Error in nexus2stype - can''t close ', filename]);
        ok=0;
    end 
end

% need to put any cognate data into tree structure
if ok
    % check that we have both tree and cognate data
    if treeread && ~isempty(content.language)
        [s,ok] = mergetreedata(s,content.array,content.language);
        if ~ok 
            disp('Error in nexus2stype: tree and data incompatible')
        end
    end
end


function [str,ok] = stype2nexus(s,comments,option,true,clades)

GlobalSwitches;

% [str,ok] = stype2nexus(s,comments,option,true,clades) has input variables
% S - a tree structure;
% Comments - a string which is included as a comment below the header in the nexus file
% Option - takes values 'DATA' 'TREE' 'BOTH' and 'TRUE' which determine
%   whether to include only the tree block, only the data block, both tree and data or
%   tree,data and synthesize block resp.
%   if option = 'TRUE' the structure True with fields mu, br, vocabsize, lambda, theta
%   must be included
% Clades - a cell array of clades.  If this is non-empty, a CLADES block will be written
%
% Returns the string str, a nexus file ready to be saved to file.
% ok == 0 if there was an error.

ok =1;
dodata = 0;
dotree = 0;
dotrue = 0;

if nargin == 5
    doclades = ~isempty(clades) & iscell(clades);
else
    doclades = 0;
end

headerstr = '';datastr = '';cladestr = '';truestr='';treestr='';
switch upper(option)
case 'BOTH'
    dodata = 1;
    dotree = 1;
case 'DATA'
    dodata = 1;
case 'TREE'
    dotree = 1;
case 'TRUE'
    if ~exist('true','var')
        disp('Error in stype2nexus: too few variables passed to write out true data')
        ok = 0;
    end
    dodata = 1;
    dotree = 1;
    dotrue = 1;
otherwise
    disp(['Error in stype2nexus: option ' option ' not recognised'])
    ok = 0;
end

if ok
    leaves = find([s.type]==LEAF);
    if isempty(leaves)
        disp('Error in stype2nexus - tree has no leaves')
        ok = 0;
    end
    % write header
    headerstr = '#NEXUS \n';
    headerstr = [headerstr '[' comments ']\n\n'];
end

if ok && dodata
% write data block
datastr = writedatablock(s,leaves);
end

if ok && doclades
    % write clades block
    cladestr = 'BEGIN CLADES;\n\n';
    for i = 1:length(clades)
        cladestr = [cladestr writeclade(clades{i})]; %#ok<AGROW>
    end
    cladestr = [cladestr 'END;\n\n'];
end

if ok && dotree
% write tree block
treestr = 'BEGIN TREES;\n\n';
treestr = [treestr 'tree treename = [&R] ' wnextree(s,find([s.type]==ROOT)) ';\n\n'];
treestr = [treestr 'END;\n\n'];
end

if ok && dotrue
% write synthesize block
truestr = 'BEGIN SYNTHESIZE;\n\nPARAMETERS ';
vals = [true.mu,true.br,true.lambda, true.theta,  true.vocabsize, true.beta ]; % LUKE added beta 04/09/2016.
truestr = [truestr sprintf('mu = %1.9e\nborrowrate = %1.9e\nlambda = %1.9e\ntheta = %1.9f\nvocabsize = %1.0f\nbeta = %1.9e;\n\n',vals)];
truestr = [truestr 'END;\n\n'];
end

str = sprintf([headerstr datastr cladestr treestr truestr]);



%----------------------------------------------------------------------
function str = writedatablock(s,leaves)

GlobalSwitches;

% write data block
str = 'BEGIN DATA;\n\n';
ntax = length(leaves);
nchar = length(s(leaves(1)).dat);
str = [str 'DIMENSIONS NTAX=' num2str(ntax) ' NCHAR=' num2str(nchar) ';\n'];
str = [str 'FORMAT MISSING=? GAP=-  INTERLEAVE ;\n\n'];
str = [str 'MATRIX \n\n'];
% get all data as a single matrix
datamatrix = reshape(sprintf('%d',[s(leaves).dat]),nchar,ntax)';
datamatrix = cellstr(datamatrix);
datamatrix = strrep(datamatrix,num2str(MIX),'?');
datamatrix = char(datamatrix);
% make list of names that can be tacked on in front of data
names = char({s(leaves).Name});
% add whitespace to end of names
names = [ names char(zeros(ntax,5)+32) ];
% make column of newline characters to tack on end of data
newline=sprintf('\n');
endcol(1:ntax,:) = newline;
% output in interleaved blocks of 100 characters
written = 0; chunksize = 100; chunks = floor(nchar/chunksize);
remain = nchar - chunksize*chunks;
if chunks==0
    % Luke 21/1/21: when nchar < chunksize, this writes a section with only taxa
    % names and no data causing the TraitLab to fail; actual data written below
    % as remain ~= 0
    % towrite = [names endcol]';
    % str = [str towrite(:)' '\n\n'];
else
    for i=1:chunks
        towrite = [names datamatrix(:,(written+1):(i*chunksize)) endcol]';
        str = [str towrite(:)' '\n\n']; %#ok<AGROW>
        if i==chunks && remain == 0
            str = str(1:end-4);
        end
        written = written + chunksize;
    end
end
if remain ~= 0
    towrite = [names datamatrix(:,(written+1):nchar) endcol]';
    str = [str towrite(:)' ';\n' ];
else
    str = [str ';\n'];
end

str = [str 'END;\n\n'];

%----------------------------------------------------------------------
function str = writeclade(c)

str = 'CLADE  ';

% write name
if ~isempty(c.name)
str = [str sprintf('NAME = %s',c.name)];
end
% write out root times
if ~isempty(c.rootrange)
    str = [str sprintf('\nROOTMIN = %1.0f ROOTMAX = %1.0f',c.rootrange)];
end
% write out adamrange
if ~isempty(c.adamrange)
    str = [str sprintf('\nORIGINATEMIN = %1.0f ORIGINATEMAX = %1.0f',c.adamrange)];
end
if isempty(c.language)
    str ='';
else
    % write out member taxa
    allnames = char(c.language);
    commas = char(44*ones(length(c.language),1));
    commas(end) = ';';
    spaces = char(32*ones(length(c.language),1));
    allnames = [allnames commas spaces]';
    str = [str sprintf('\nTAXA = %s\n\n',allnames)];
end

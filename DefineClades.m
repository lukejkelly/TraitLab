function clade=DefineClades(activeclades)

% function clade=DefineClades(activeclades)
%
% Indo-European clades of interest
% GKN 26/6/03
% 
% nargin=1
% if activeclades cell array of strings is present
% data is restricted to clades in cell array (so data
% from clades missing from activeclades is not imposed


clade=DefineCladesData;

if nargin>0
    nc=size(clade,2);
    for t=1:nc
        if ~any(strcmp(activeclades,clade{t}.name))
            clade{t}.rootrange=[];
            clade{t}.adamrange=[];
        end
    end
end
function [nodeTimes] = getConstraintNodeTimes(pathToTrees, constraintCA, ...
                                              sampInds, refNodeName)
% Return node times of most recent common ancestor of nodes
%   * pathToTrees:  path to output file
%   * constraintCA: cell array of leaf node names; e.g. {'1'} or {'1', '2'}
%   * sampInds:     indices of samples we want returned
%   * refNodeName:  reference node whose time we compare against; this can be
%                   empty if constraintCA has more than one entry
% We need to specify a reference node as when the tree is written to file, the
% node times times are shifted so that the most recent node has time 0.

% Global variables.
GlobalSwitches; GlobalValues;

% Reading trees.
trees = readalltrees(pathToTrees);

% Empty vector of node times.
nodeTimes = zeros(length(sampInds), 1);

% Cycling through trees.
for i = 1:length(sampInds)

  state = tree2state(rnextree(trees{sampInds(i)}));
  node  = findCA(state, constraintCA);
  c     = state.tree(node).time;

  % Reference node time
  if length(constraintCA) == 1
    d = state.tree(state.leaves(cellfun(@(X) strcmp(X, refNodeName), ...
                                        {state.tree(state.leaves).Name}))).time;
  else
    d = 0;
  end

  nodeTimes(i) = c - d;
end

end

% Find MRCA in state of set of leaf nodes CA.
function [caNode] = findCA(state, CA)

% Number of leave nodes to find ancestor of.
lCA = length(CA);

if lCA == 1
  for j = state.leaves
    if strcmp(state.tree(j).Name, CA{1})
      caNode = j;
    end
  end
else
  % Getting indices of CA nodes in tree.
  childInds = zeros(size(CA));
  for j = 1:lCA
    for k = state.leaves
      if strcmp(state.tree(k).Name, CA{j})
        childInds(j) = k;
      end
    end
  end

  % Finding MRCA.
  ancInds = unique([state.tree(childInds).parent]);
  while ~all(ancInds == ancInds(mod((1:length(ancInds)), length(ancInds)) + 1))
    inds = find([state.tree(ancInds).time] == min([state.tree(ancInds).time]));
    ancInds(inds) = [state.tree(ancInds(inds)).parent];
    ancInds       = unique(ancInds);
  end
  caNode = ancInds;
end

end

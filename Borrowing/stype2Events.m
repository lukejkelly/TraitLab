function [ tEvents, rl ] = stype2Events( s )
% Takes as input the struct s detailing the nodes. Puts a right-to-left
% order on the tree topology then creates a struct tEvents to be used to
% find the likelihood.
%   - tEvents = { L, time, type, location, K }.
%   - event types =
%       - 1, branching event.
%       - 2, catastrophe.
%       - 3, branch death.
%   - eventLocation = branch index (counting from right to left).
% Catastrophes on the Adam-root branch give spurious output.

% Declaring global variables
global LEAF ANST ROOT

% Correcting slight discrepancies in leaf times.
for k = find( [s.type] == LEAF & [s.time] < min([s.time]) + 0.01 )
    
    s(k).time = min([s.time]);
    
end

% Node times, starting with the leaves and finishing with root.
times = unique( [ s.time ] );

% List of visited nodes.
visited = zeros( size(s) );

% Ordered list of branches in tree at each time. The cell rtol{k} gives
% the right-left ordering of the branches/nodes at times(k). 
rtol = cell( size(times) );

% Starting location.
loc = find( [s.type] == LEAF, 1 );

% We now visit each node along a path, starting at node loc. loc is
% considered to be the right-most node at the time at which it occurs and
% the other nodes are added in the order visited and relative to loc. This
% induces a topologicial representation of the tree, which in turn induces
% the patterns and the basis for solving the differential equations.

while ~all(visited) % While each node has yet to be visited.
    
    if visited(loc) == 0 % If node hasn't been visited, record its details.
    
        % Add location to list of visited nodes.
        visited(loc) = loc;

        % Time slots on which edge between node and its parent existed.
        if s(loc).type == LEAF

            % For branches which die, we record their order as if they were
            % still evolving so that we can treat the patterns properly.
            slots = find( times < s( s(loc).parent ).time );
                        
        elseif any( s(loc).type == [ANST, ROOT] )

            slots = find(times >= s(loc).time & times < s(s(loc).parent).time);

        else % Adam node.
            
            slots = find( times >= s(loc).time );

        end
        
        % Add location to branch order for each time slot.
        for k = slots
            
            rtol{k} = [ rtol{k}, loc ];
            
        end

    else % If node has already been visited, move to a neighbouring node.
        % Is this 'else' statement necessary?
        
        % If at a leaf, move upwards.
        if s(loc).type == LEAF

            loc = s(loc).parent;

        % If at a node whose children have both been visited, move up.
        elseif all( visited(s(loc).child) ) % all( ismember( s(loc).child, visited ) )

            loc = s(loc).parent;

        % If neither, move down.
        else % Move to second child if first has already been visited.

            if visited(s(loc).child(1)) % ismember( s(loc).child(1), visited )

                loc = s(loc).child(2);

            else

                loc = s(loc).child(1);

            end

        end % End of movement.
    
    end % End of action.
    
end % End of while loop.

% Right-to-left ordering at leaves, used to calculate likelihood.
rl = rtol{1};

% Now we add branch information for each node - namely, the branch on which
% a node arises / disappears / etc.
for k = 1:length(s)
    
    s(k).branch = find( rtol{ s(k).time == times } == k );
    
end

% Lists of branches which are dead at each event time.
K = cell( length(times), 1 );

% We populate the lists of dead branches for each topological event time,
% that is branching events or branch deaths.
for k = find( [s.type] == LEAF & [s.time] > min([s.time]) )
    
    % For the time at which a branch dies, we record its index and add it to K.
    ind = find(times == s(k).time);
    K{ind} = [ K{ind}, find( rtol{ind} == k ) ];
    
    % For the event times thereafter, we record the branch index just after
    % the event (the rtol order just before the next event). We don't
    % require the list of dead branches at the leaves as this is not the
    % start of an interval on which we solve the differential equations.
    slots = find( times < s(k).time & times > min([s.time]) );
    
    for j = slots
        
        K{j} = [ K{j}, find( rtol{j - 1} == k ) ];
        
    end    
    
end

% Having obtained the branch ordering of the nodes, we now need the same
% for the catastrophes and the times at which they occur.
% Create an array of catastrophe information. Three columns:
%   - Time,
%   - Index of node such that catastrophe occurs on branch linking node to
%         its parent,
%   - Catastrophe branch order in right-left ordering.

% Indices of nodes with catastrophes on thir parent branches.
cat_inds = find( cellfun( @any, {s.catloc} ) );

% Number of catastrophes on the tree.
n_cat = 0;

for k = cat_inds % This needs to change.
    
    n_cat = n_cat + length( s(k).catloc );
    
end

% Empty array of catastrophe information.
cats = zeros( n_cat, 3 );

% Cycle through nodes with catastrophes on the branches linking them to
% their parents.
i = 1;

for k = cat_inds
    
    % Number of catastrophes on the node.
    len = length( s(k).catloc );
    
    % Populating time and branch index columns of cats.
    cats( i:(i + len - 1), 1 ) = s(k).time + ...
        s(k).catloc * ( s( s(k).parent ).time - s(k).time );
    cats( i:(i + len - 1), 2 ) = k;
    
    % Populating right-left branch order column of cats. If a catastrophe
    % occurs on the branch linking node i to its parent, <i, pa(i)>, then
    % the right-left order is given by that of i's order at the next
    % topological event on the tree.
    
    j = 1;
    
    while j <= len

        time_j = max( times( times < cats(i + j - 1, 1) ) );
        
        cats(i + j - 1, 3) = find( rtol{ times == time_j } == ...
            cats(i + j - 1, 2) );

        j = j + 1;
        
    end

    i = i + len;

end

% Creating treeEvents array with three columns:
%   - Time,
%   - Type,
%   - Location.

% Some variables to be used in indexing.
eventTime = 1;
eventType = 2;
eventLocation = 3;

% Numbers of offset leaves, branching events and catastrophes.
n_uleaves = sum( [s.type] == LEAF & [s.time] > min(times) );
n_nodes = sum( [s.type] == ANST ) + 1; % +1 for root node.

% Empty matrix of tree events.
treeEvents = zeros( n_nodes + n_cat + n_uleaves + 1, 3 );

% Counter
i = 1;

% Adding entry for leaf times to matrix.
treeEvents(i, eventTime) = min(times);
treeEvents(i, eventType) = 0;
treeEvents(i, eventLocation) = 0;

i = i + 1;

% Adding branching events to matrix.
for k = find( [s.type] == ANST | [s.type] == ROOT )
    
    treeEvents(i, eventTime) = s(k).time;
    treeEvents(i, eventType) = 1;
    treeEvents(i, eventLocation) = s(k).branch;

    i = i + 1;
    
end

% Adding catastrophes to matrix.
treeEvents(i:(i + n_cat - 1), eventTime) = cats(:, 1);
treeEvents(i:(i + n_cat - 1), eventType) = 2;
treeEvents(i:(i + n_cat - 1), eventLocation) = cats(:, 3);

i = i + n_cat;

% Adding branch deaths.
for k = find( [s.type] == LEAF & [s.time] > min(times) )

    treeEvents(i, eventTime) = s(k).time;
    treeEvents(i, eventType) = 3;
    treeEvents(i, eventLocation) = s(k).branch;
    
    i = i + 1;
    
end

% Sort by eventTime, starting with the first branching event.
treeEvents = sortrows(treeEvents, - eventTime);

% Adding a column for the number of lineages during each interval.
treeEvents = [ 1 + cumsum( treeEvents(:, eventType) == 1 ), treeEvents ];

% Creating a struct, tEvents to return the above array along with the
% details of the dead leaves and evolving branches.
tEvents = struct( ...
    'L', {[]}, ...
    'time', {[]}, ...
    'type', {[]}, ...
    'loc', {[]}, ...
    'K', {[]} );

[len, ~] = size(treeEvents);

for k = 1:len
    
    tEvents(k).L = treeEvents(k, 1);
    tEvents(k).time = treeEvents(k, 2);
    tEvents(k).type = treeEvents(k, 3);
    tEvents(k).loc = treeEvents(k, 4);
    
    % Add the list of dead branches at the time of the event (if a branch
    % birth or death) or at the preceding such event (if a catastrophe).
    tEvents(k).K = K{ find(times >= tEvents(k).time, 1) };

end

end
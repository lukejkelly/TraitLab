function [ patternMeans, patternVars ] = generateTestPatterns( tree, rates, nRepeats )
%Generates data on a phylogenetic tree, without borrowing
%   tree is the tree structure, arranged as a matrix.
%      size (number of vertices, 5), where each row takes the form:
%      [VertexID ParentID LeftChildID RightChildID Time]
%   rates is a struct containing the birth death and borrowing parameters
%   repeats is how many times to do the process to build up distributions
%   patternMeans are the observations

if rates.borrowing ~= 0
    fprintf('**Warning: these simulations do not support borrowing.\n')
end

%Define some constants
setParameters;

%Plot tree to get things started
%plotTree(tree, 'labels', true);

%Record tree structure
LtoR = parseTreeStructure(tree);
nVertices = size(tree, 1);
L = nVertices + 1;

%No. of possible site patterns
nPatterns = 2^L - 1;
patternMeans = zeros(nPatterns, nRepeats);

%Sorting system
ind = 1:nPatterns;  %All site pattern numbers
%Sort site patterns into sets.
binSorted = sortPatterns( de2bi(ind) );

%Loop through tree, from oldest first solving equations
%Sort the tree by time, oldest first
[~,order] = sortrows(tree, -TIME);

%Repeat simulations to build up distributions
for i = 1:nRepeats
    
    traits = { };
    
    %Generate initial set of traits at oldest vertex
    traits{ order(1) } = 1:poissrnd(rates.birth / rates.death);
    maxTrait = max( traits{order(1)} );
    if isempty(maxTrait)
        maxTrait = 0;
    end
    
    %Loop through the rest of the vertices, generating traits
    for j = 2:nVertices
        
        vertexRow = order(j);
        vertex = tree(vertexRow, ID);
        
        parent = tree(vertexRow, PARENT);
        parentRow = find(tree(:, ID) == parent);
        
        %Calculate time to next vertex
        tParent = tree(parentRow, TIME);
        tVertex = tree(vertexRow, TIME);
        t = tParent - tVertex;
        
        %Generate new traits
        [traits{vertexRow}, maxTrait] = updateTraits(traits{parentRow}, ...
            rates, t, maxTrait);
        
    end
    
    %Then loop through leaves, generating traits
    leafTraits = { }; l = 1;
    for j = 1:nVertices
        
        %Left to right order
        vertex = LtoR{end}(j);
        vertexRow = find(tree(:, ID) == vertex);
        t = tree(vertexRow, TIME);
        
        if tree(vertexRow, LEFTCHILD) == LEAF
            
            %Update traits
            [leafTraits{l}, maxTrait] = updateTraits(traits{vertexRow}, ...
                rates, t, maxTrait);
            l = l+1;
            
        end
        
        if tree(vertexRow, RIGHTCHILD) == LEAF
            
            %Update traits
            [leafTraits{l}, maxTrait] = updateTraits(traits{vertexRow}, ...
                rates, t, maxTrait);
            l = l+1;
            
        end
        
    end
    
    %Finally turn into site patterns and store
    patterns = zeros(maxTrait, L);
    for j = 1:L
        patterns( leafTraits{j}, j ) = 1;
    end
    
    %Loop through site patterns
    for j = 1:nPatterns
        
        %Compare against observed trait
        for k = 1:maxTrait
            
            %If the same, store
            if isequal( binSorted(j,1:end-1), patterns(k,:) )
                patternMeans(j, i) = patternMeans(j, i) + 1;
            end
            
        end
        
    end
    
end

patternVars = var(patternMeans, 0, 2);
patternMeans = mean(patternMeans, 2);

end

function [ newTraits, maxTrait ] = updateTraits( oldTraits, rates, t, maxTrait )
%Given a set of traits, rates and a time of existence calculates a new set
%of traits. First culls existing traits, then generates new ones

%Copy traits from parent set to here
survivingTraits = oldTraits;

%Work out probability of survival of traits to here
pSurvive = exp(-rates.death * t);

%Cull traits from parent set
survivingTraits( rand(size(survivingTraits)) > pSurvive ) = [ ];

%Add new traits
mBorn = (rates.birth / rates.death) * (1 - exp(-rates.death*t));
bornTraits = (1:poissrnd(mBorn)) + maxTrait;
%Update maximum existing trait
maxTrait = max([maxTrait bornTraits]);
if isempty(maxTrait)
    maxTrait = 0;
end

%Save traits
newTraits = [survivingTraits bornTraits];

end


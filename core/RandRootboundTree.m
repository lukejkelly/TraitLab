function s = RandRootboundTree(clade,langname, rootmax)
% builds a random tree using the coalescent tree building model
% so that all clades are observed and the times are scaled so that 
% the age constrains on the clades are observed
% assumes that if clade a and b have anything in common, a contains b or b contains a 

% turn clades into a structure array
x=[clade{:}];

% make sure that there is a clade containing all languages
biggestclade = 1;
nclade = length(x);
rootmin = 0;
for i = 1:nclade
    if length(x(i).language) > length(x(biggestclade).language)
        biggestclade = i;
    end
    if x(i).rootrange(1)>rootmin
        rootmin = x(i).rootrange(1);
    end
end

if length(x(biggestclade).language) < length(langname)
    % need to create the super clade
    nclade = nclade +1;
    % make sure vector is a row vector
    if size(langname,1) > size(langname,2)
        x(nclade).language = langname';
    else
        x(nclade).language = langname;
    end
    x(nclade).name = sprintf('Clade%f',nclade);
    x(nclade).rootrange = [rootmin rootmax];
    x(nclade).adamrange = [];
else
    % move super clade to the end
    x(end+1) = x(biggestclade);
    x(biggestclade) = [];
end

% figure out a building order for the clades
contains = zeros(nclade);
for i = 1:nclade
    for j = i+1:nclade
        incommon = intersect(x(i).language,x(j).language);
        if ~isempty(incommon)
            if length(x(i).language)<length(x(j).language)
                contains(j,i) = 1;
            else
                contains(i,j) = 1;
            end
        end
    end
end
ncontain = sum(contains,2)';

% make sure clade is only listed as a subclade of one other
nsup = sum(contains,1);
for i = 1:nclade
    if nsup(i)>1
        supi = find(contains(:,i));
        smallest = supi(1);
        for k = 1:nsup(i)
            if length(x(supi(k)).language) < length(x(smallest).language)
                smallest = supi(k);
            end
        end
        contains(:,i) = zeros(nclade,1);
        contains(smallest,i) = 1;
    end
end

% if a clade contains no others, we can build it now
% if it does contain others, we must build them first

cladetree.s = [];
[cladetree(1:nclade).roottime] = deal(0);

done = zeros(1,nclade);

% build clade trees in order of how many subclades each has starting from zero
for i = 0:max(ncontain)
    % make those with i subclades
    j = find(ncontain == i);
    for k = j
        % make up a list of languages contained with in this clade but not any subclades
        n = find(contains(k,:));
        lang = x(k).language;
        for m = n
            lang = setdiff(lang,x(m).language);
        end
        time = zeros(1,length(lang));
        for m = n
            lang(end+1) = {x(m).name};
            time(end+1) = cladetree(m).roottime;
        end
        [cladetree(k).s cladetree(k).roottime] = BuildOffsetTree(lang,time,x(k).rootrange(1),x(k).rootrange(2));
        for m = n
            ileaf = find(strcmp({cladetree(k).s.Name},x(m).name));
            cladetree(k).s = jointree(cladetree(k).s,cladetree(m).s,ileaf);
        end
    end
end

s = cladetree(end).s;
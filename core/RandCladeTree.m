function s = RandCladeTree(clade,langname, rootmax,rho);
% function s = RandCladeTree(clade,langname, rootmax);
% builds a random tree using the coalescent tree building model
% so that all clades are observed and the times are scaled so that 
% the age constrains on the clades are observed
% assumes that if clade a and b have anything in common, a contains b or b contains a 

%keyboard;
%%
% turn clades into a structure array
x=[clade{:}];
nclade = length(x);

% make sure that there is a clade containing all languages
% check which clades need imposition of rootmax bounds
% check which clades have bounds on adam node and need to be converted to root bounds
biggestclade = 1;
adambound = zeros(1,nclade);
noroot = 0;
for i = 1:nclade
    if length(x(i).language) > length(x(biggestclade).language)
        biggestclade = i;
    end
    if ~isempty(x(i).rootrange)
        if isinf(x(i).rootrange(2))
            x(i).rootrange(2) = rootmax;
        end
        noroot = 0;
    else
       noroot = 1; 
    end
    if ~isempty(x(i).adamrange)
        adambound(i) = 1;
        if isinf(x(i).adamrange(2))
            x(i).adamrange(2) = rootmax; 
        end
    else
        if noroot
            x(i).rootrange = [0 rootmax];
        end
    end
end

%%

% build a clade with all languages if it doesnt already exist
if length(x(biggestclade).language) < length(langname)
    % need to create a clade containing al languages
    nclade = nclade +1;
    x(nclade) = pop('clade');
    % make sure vector is a row vector
    if size(langname,1) > size(langname,2)
        x(nclade).language = langname';
    else
        x(nclade).language = langname;
    end
    x(nclade).name = sprintf('Clade%f',nclade);
    x(nclade).rootrange = [0 rootmax];
    biggestclade = nclade;
    adambound(nclade) = 0;
    nobound(nclade) = 0;
else
    % if maximal clade has bound on adam, change it to a bound on the root
    % TODO - not sure why we replace the rootrange with adamrange - 
    %      - if user gave both and are unequal unclear why adamrange is
    %      - used. GKN, 18/8/10
    if adambound(biggestclade)
        x(biggestclade).rootrange = x(biggestclade).adamrange;
        x(biggestclade).adamrange = [];
        adambound(biggestclade) = 0;
    end
end
%%
% figure out a building order for the clades by finding the 
% partial ordering of clades
% make sure clade is only listed as a subclade of one other
[contains,ncontain] = cladeorder(x);

%%
% check that all bounds are consistent - ie bounds on a clade are >= than bounds on its subclade
for i = max(ncontain):-1:0
    % check those with i subclades
    j = find(ncontain == i);
    for k = j
        %TODO - this looks wrong - the upper bound on k's root range
        %     - is below the upper bound on k's adam, and the lower
        %     - upper bound bounds the adam and root of any clade beneath
        %     - so I suggest
        %upperbound = x(k).rootrange(2); % irrespective of adamrange
        if adambound(k)
            upperbound = x(k).adamrange(2);
        else
            upperbound = x(k).rootrange(2);
        end
        n = find(contains(k,:));
        for m = n
            if adambound(m)
                if x(m).adamrange(2) > upperbound
                    x(m).adamrange(2) = upperbound;
                end
            else
                if x(m).rootrange(2) > upperbound
                    x(m).rootrange(2) = upperbound;
                end
            end
        end
    end
end

%%
% go through clades with adambounds creating a super clade for each of them
% which includes one language or clade from existing superclade
% and has a bound on the root
for i = max(ncontain):-1:0
    % check those with i subclades
    j = find(ncontain == i);
    for k = j
        if adambound(k) 
            % need to make superclade 
            % get clade which is smallest clade containing this one
            nextup = find(contains(:,k));
            % make up a list of languages contained in parent clade but not any subclades
            n = find(contains(nextup,:));
            lang = x(nextup).language;
            for m = n
                lang = setdiff(lang,x(m).language);
            end
            % if lang is empty, need to choose a clade to latch on to
            if ~isempty(lang)
                % randomly choose one of the remaining languages,
                % chuck it and langs from subclade to make a new clade
                xnew = pop('clade');
                xnew.language = [x(k).language lang(ceil(rand*length(lang)))];
                xnew.rootrange = x(k).adamrange;
                xnew.name = ['Supercladeofcladecalled' x(k).name];
                % check that this clade is not just a copy of the parent
                if (length(xnew.language) == length(x(nextup).language)) & ~adambound(nextup)
                    % it is, see if we can reconcile root bounds
                    xnew.rootrange = [max([xnew.rootrange(1) x(nextup).rootrange(1)]) min([xnew.rootrange(2) x(nextup).rootrange(2)])];
                    if xnew.rootrange(1) >=  xnew.rootrange(2)
                        % protest just stop and ask that file be changed
                        error(sprintf('Cant build clade tree with given bounds due to apparent redundancy between clades %s and %s',x(k).name,x(nextup).name)); 
                    else
                        % these new bounds are ok - modify parent
                        x(nextup).rootrange = xnew.rootrange;
                        %GKN 5/1/10 next line
                        x(k).rootrange = xnew.rootrange;
                        x(k).adamrange = [];    
                        
                    end    
                else
                    % its not, add to end
                    x(end+1) = xnew;
                    adambound(end+1)=0; %RJR 17/07/08
                    %change original clade %GKN 5/1/10
                    if isempty(x(k).rootrange), x(k).rootrange = [0, xnew.rootrange(2)]; end
                    x(k).adamrange = [];   
                    
                end    
            else
                % need to choose another clade to buddy up with
                % go through clades, looking for one with either a root time less than this adam time 
                % or an adam time overlapping this one
                n(find(n==k)) = [];
                buddy = 0;
                for m = n
                    % see if bounds are compatible
                    if adambound(m)
                        range = [max([x(m).adamrange(1) x(k).adamrange(1)]),min([x(m).adamrange(2) x(k).adamrange(2)]) ];
                    else
                        range = [x(m).rootrange(1) x(k).adamrange(2)];
                    end 
                    if range(1) < range(2)
                        xnew = pop('clade');
                        xnew.language = [x([k,m]).language];
                        xnew.name = ['Supercladeofcladecalled' x(k).name 'and' x(m).name];
                        xnew.rootrange = range;
                        % check to see if this is not just a copy of the parent
                        if length(xnew.language) < length(x(nextup).language) 
                            buddy = 1;
                            if adambound(m)
                                %GKN 5/1/10
                                if isempty(x(k).rootrange), x(k).rootrange = [0, xnew.rootrange(2)]; end
                                x(m).rootrange = [0 xnew.rootrange(2)];
                            else
                                %GKN 5/1/10
                                if isempty(x(k).rootrange), x(k).rootrange = [0, rootmax]; end
                                x(m).rootrange(2) = min([x(m).rootrange(2),xnew.rootrange(2)]);
                            end
                            x(end+1) = xnew;
                            adambound(end+1)=0; %RJR 17/07/08
                            break   
                        else
                            % it is a copy of parent - check that bounds are compatible
                            xnew.rootrange = [max([xnew.rootrange(1) x(nextup).rootrange(1)]) min([xnew.rootrange(2) x(nextup).rootrange(2)])];
                            if xnew.rootrange(1) >=  xnew.rootrange(2)
                                % protest just stop and ask that file be changed
                                error(sprintf('Cant build clade tree with given bounds due to apparent redundancy between clades %s and %s',x(k).name,x(nextup).name)); 
                            else
                                % these new bounds are ok - modify parent
                                x(nextup).rootrange = xnew.rootrange;
                                % modify original adambound clade
                                %GKN 5/1/10
                                if isempty(x(k).rootrange), x(k).rootrange = [0, xnew.rootrange(2)]; end
                                x(k).adamrange = [];  
                                %modify buddy clade
                                x(m).rootrange = [0 xnew.rootrange(2)];
                                x(m).adamrange = [];  
                                buddy = 1;
                            end    
                        end
                    end
                end
                if ~buddy 
                    % no buddy - error
                    error(sprintf('Cant build clade tree with given bounds due to apparent redundancy between clades %s and %s',x(k).name,x(nextup).name))
                else
                    adambound([k,m])=0;
                end
            end
            [contains] = cladeorder(x);
        end
    end
end

%%

% figure out a building order for the clades by finding the partial ordering of clades
% make sure clade is only listed as a subclade of one other
[contains,ncontain] = cladeorder(x);

%%
% initialize tree structure - have one tree for each clade
cladetree.s = [];
[cladetree(1:nclade).roottime] = deal(0);


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
        [cladetree(k).s cladetree(k).roottime] = BuildOffsetTree(lang,time,x(k).rootrange(1),x(k).rootrange(2),rho);
        for m = n
            ileaf = find(strcmp({cladetree(k).s.Name},x(m).name));
            cladetree(k).s = jointree(cladetree(k).s,cladetree(m).s,ileaf);
        end
    end
end
%%
s = cladetree(biggestclade).s;

if length(s)~=2*length(langname)
    disp('unequal # taxa in RandCladeTree');keyboard;pause;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper functions
function [contains,ncontain] = cladeorder(x)

nclade = length(x);

% figure out a building order for the clades by finding the partial ordering of clades
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

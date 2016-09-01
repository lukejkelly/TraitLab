function [rseq,langs,cogs,NS,L,L_vals,s,catnum]=PolyMorph(br,mu,p,lambda_total,s,min_wpm,N_Meaning_Classes,local,maxdist,branchvar,classvar,simcats,rho,kappa,knowcats)
%function [rseq,langs,cogs,NS,L,L_vals,mseq]=PolyMorph(br,mu,p,lambda_total,s,min_wpm,N_Meaning_Classes,local,maxdist,branchvar,classvar)
% br is borrowing rate
% mu is death rate per word
% lambda_total is total per language birth rate
% s is tree
% min_wpm is minimum number of words per meaning
% local = ON for local borrowing = OFF for wide borrowing
% maxdist - if langs have common ancestor within maxdist then they can borrow with each other
% if rhbranchvar > 0 then rate heterogeneity across branches according to
% LogNormal dist with mean=mu, var=branchvar
% if rhclassvar > 0 then rate heterogeneity across obs classes according to
% LogNormal dist with mean=mu_i, var=classvar where mu_i is rate on
% branch (just mu if branchvar == 0)

GlobalSwitches

% get leaves
leaves=find([s.type]==LEAF);
% count leaves
NS=length(leaves);
% get root index
Root=find([s.type]==ROOT);

if ~exist('N_Meaning_Classes','var')
    if min_wpm==1
        N_Meaning_Classes=input('N_Meaning_Classes ');
        %N_Meaning_Classes=round(lambda_total/mu);
    else
        N_Meaning_Classes=1;
    end
end
if ~exist('local','var')
    local = OFF;
end
if local==ON & ~exist('maxdist','var')
    error('Local borrowing called for but maxdist not defined');
end
if local==ON & maxdist<0
    error('maxdist must be >= 0 for local borrowing');
end
if local == OFF
    % widespread borrowing - set maxdist > tree height
    maxdist = s(Root).time +1;
end
if ~exist('branchvar','var')
    branchvar = 0;
    classvar = 0;
end
if ~exist('classvar','var')
    classvar = 0;
end
if ~exist('simcats','var')
    simcats=0;
end
if ~exist('knowcats','var')
    knowcats=0;
end


% divide total birth rate evenly between the meaning classes
lambda=lambda_total/N_Meaning_Classes;

% L_total total number of words regardless of meaning class
% L_vals are cut-off points for words in different meaning classes
L_total=0;L_vals=[];

% see if there is rate heterogeneity across branches
if branchvar == 0
    % store death rate for branch in its child vertex
    [s.mub] = deal(mu);
elseif branchvar < 0
    error('branchvar (rate variance across branches) is negative')
else
    % have rate het on branches
    % generate rate for each branch
    allmu = randln(mu,mu*branchvar,1,length(s));
    % store rate for branch in its child vertex
    for i=1:length(s)
        s(i).mub = allmu(i);
  %      disp(sprintf('On edge %g, mu(edge)/mu(mean) = %g ',i,allmu(i)/mu));
    end
end

%generate catastophes
if simcats 
    if ~knowcats
        catnum=zeros(1,2*NS);
        catlist=[];
        for i=1:length(s)
            s(i).cat=[];
            if s(i).type<=ANST
                ncat(i)=MyPoisson(rho*(s(s(i).parent).time-s(i).time));
                catnum(i)=ncat(i);
                for j=1:ncat(i)
                    s(i).cat(j)=s(i).time+rand*(s(s(i).parent).time-s(i).time);
                    catlist(:,end+1)=[i;s(i).cat(j)];
                end
            end
        end

        tc=-log(1-kappa)/mu;
        disp(['Simulated ' num2str(sum(ncat)) ' catastrophes.']);
    else % old catastrophes are already in the tree; use those
        catnum=zeros(1,2*NS);
        catlist=[];
        for i=1:length(s)
            if s(i).type<=ANST
                catnum(i)=length(s(i).cat);
                for j=1:catnum(i)
                    catlist(:,end+1)=[i;s(i).cat(j)];
                end
            end
        end
        tc=-log(1-kappa)/mu;
        ncat=length(catlist);
    end
    if ncat==0, simcats=0; end
else
    catnum=[];
end
            


for m=1:N_Meaning_Classes
    if N_Meaning_Classes>1 && mod(m,20) == 0
        disp(sprintf('Simulations running: up to observation class %d of %d',m,N_Meaning_Classes));
    end
    
    if classvar == 0
        mucl = 1;
    else
        mucl = randln(1,classvar,1);
  %      disp(sprintf('In observation class %g, mu(obs class)/mu(mean) = %g ',m,mucl));
    end

    % sample equilibrium number of words at root of tree
    L=MyPoisson(lambda/(s(Root).mub*mucl));
    check=0;
    while L<min_wpm
        % need to resample until there are at least min_wpm at the root
        L=MyPoisson(lambda/(s(Root).mub*mucl)); check=check+1;
        if check>50,
            % poisson sampler inefficient - just assign a single word at the root
            disp(sprintf('L-generation inefficient in PolyMorph() - simulating class size zero.\n'));
            disp(sprintf('Assigning L=1 WARNING - continue if you know what you are doing.\n'));
            L=1;
            keyboard;
        end
    end
    % assign words to root (words names are just numbers so first L words are 1:L)
    seq(Root).words=1:L;
    % keep count of words at root node
    seq(Root).L=L;
    if seq(Root).L~=numel(seq(Root).words), keyboard, end

    % cnt is current time
    cnt=s(Root).time;
    % n is list of languages during current time slice (lang names are node at bottom of edge)
    n=[];
    % copy words from root down to child accounting for effects to cladegenic loss
    for c=s(Root).child
        if p < 1
            % have cladagenic loss
            if min_wpm==0
                [seq(c).words,seq(c).L]=BranchVocab(seq(Root).words,p);
            elseif min_wpm>=1
                [seq(c).words,seq(c).L]=BranchVocabPolyMorph(seq(Root).words,p,min_wpm);
            else
                error('Weirdness in PolyMorph(): min_wpm should be a whole number (for the moment, 0 or 1).');
            end
        else
            % no cladagenic loss - just directly copy parent language
            seq(c).words = seq(Root).words;
            seq(c).L = seq(Root).L;
             if seq(c).L~=numel(seq(c).words), keyboard, end
        end
    end
    % add child nodes to list of current langs
    newn=[s(Root).child];
    n=[n,newn];
    % nb is list of langs that are within borrowing distance of another
    nb = n;
    % next time is first branching event among current langs
    nnt=max([s(n).time]);
    
    %next catastrophe time
    if simcats
        tempcatlist=catlist; %this list will be shortened as catastrophes are encountered
        [nc idx]=max(tempcatlist(2,:));
        nci=tempcatlist(1,idx);
        
    else
        nc=-Inf;
        foundcat=0;
    end

    % dist keeps track of distance between language - infinite if
    dist = repmat(inf,2*NS-1);
    dist(n(1),n(2)) = 0;
    dist(n(2),n(1)) = 0;

    while 1
        % nw is total number of words  (same word in 2 different langs counts as 2)
        nw=sum([seq(n).L]);
        % num words for borrowing is all words within maxdist of another language
        nwb = sum([seq(nb).L]);

        if local ==OFF & nw ~=nwb
            error('Local borrowing must be on as nw ~= nwb but local == OFF.')
        end
        numlang=length(n);


        % T is length of interval where we have this set of languages
        T=cnt-nnt;
        % y is amount of time already spent in this interval
        y=0;
        % tdist is next time that some languages become too distant to borrow from each other
        tdist = min(min(maxdist - dist(isfinite(dist))));

        while 1
            % total rate that events happen is sum(localdeathrate*nwordsinlang)+nborrowwords*borrowrate + nlang*birthrate
            rate = mucl * sum([seq(n).L] .* [s(n).mub]) + nwb*br + numlang*lambda;
            % draw exponential rand at total rate
            yinc = -log(rand)/rate;
            y = y+yinc;
            
            % check whether we have gone through a catastrophe
            if simcats, foundcat=((y>cnt-nc) && (nc>nnt)); end
            while foundcat
                % deaths at lang nci
                %was rand(seq(nci).L,1)<1-exp(-s(nci).mub*tc) but that is
                %wrong, I think it should be as below - GKN 1/4/11
                nseq= seq(nci).words(rand(seq(nci).L,1)<exp(-s(nci).mub*tc));
                                
                % births at lang nci
                ncb=MyPoisson(kappa*lambda/mu);
                nseq=[nseq, L+1:L+ncb];
                L=L+ncb;
                
                if length(nseq)<min_wpm % keep one at random
                    nseq=seq(nci).words(ceil(seq(nci).L*rand));
                end
                
                nw=nw-seq(nci).L+length(nseq);
                
                if ~isempty(nb) && any(nb==nci)
                    nwb = nwb-seq(nci).L+length(nseq);
                end
                
                seq(nci).L=length(nseq);
                seq(nci).words=nseq;
                 if seq(nci).L~=numel(seq(nci).words), keyboard, end
                
                tempcatlist(:,idx)=[];
                
                if ~isempty(tempcatlist)
                    [nc idx]=max(tempcatlist(2,:));
                    nci=tempcatlist(1,idx);
                    foundcat=((y>cnt-nc) && (nc>nnt));
                else
                    nc=-Inf;
                    foundcat=0; 
                end
            end
            
            
            % if gone beyond interval, stop
            if y>T, dist = dist - y + yinc + T; break; end;

            if y > tdist
                % gone beyond tdist, adjust dist, nwb, y and try again
                dist = dist - y + yinc + tdist;
                y = tdist;
                dist(dist == maxdist) = inf;
                tdist = y+min(maxdist - dist(isfinite(dist)));
                nb = near(dist, maxdist);
                nwb = sum([seq(nb).L]);
            else
                dist = dist + yinc;

                % chose another rand to see whether event is borrowing, birth or death
                r=rand;
                if r<nwb*br/rate
                    % its a borrowing
                    % choose language which lends word weighted by number of words in the language
                    lngout=nb(disample([seq(nb).L]./nwb));
                    % choose a word from that language
                    wd=ceil(rand*seq(lngout).L);
                    % choose languge which borrows word
                    blang = find(dist(lngout,:)<maxdist);
                    lngin=blang(ceil(rand*length(blang)));
                    % add if not already present
                    if ~any(find(seq(lngin).words==seq(lngout).words(wd)))
                        % word not in language
                        % add it
                        seq(lngin).words = [seq(lngin).words,seq(lngout).words(wd)];
                        % update nw, seq.L, nwb
                        seq(lngin).L=seq(lngin).L+1;
                         if seq(lngin).L~=numel(seq(lngin).words), keyboard, end
                        nw = nw+1;
                        nwb = nwb+1;
                    end
                elseif r<(nwb*br + numlang*lambda)/rate
                    % its a birth
                    % choose lang UAR
                    lng=n(ceil(rand*numlang));
                    % pop new word of stack
                    L=L+1;
                    seq(lng).words=[seq(lng).words,L];
                    seq(lng).L=seq(lng).L+1;
                     if seq(lng).L~=numel(seq(lng).words), keyboard, end
                    nw=nw+1;
                    if ~isempty(nb) && any(nb==lng)
                        nwb = nwb+1;
                    end
                else
                    % it's a death
                    % choose language weighted by number of words and local
                    % death rate
                    wght = [seq(n).L] .* [s(n).mub];
                    lng=n(disample(wght/sum(wght)));
                    % check that we are not going below min words per meaning
                    if seq(lng).L>min_wpm
                        % choose a word UAR
                        wd=ceil(rand*seq(lng).L);
                        % get rid of it
                        seq(lng).words(wd)=[];
                        seq(lng).L=seq(lng).L-1;
                         if seq(lng).L~=numel(seq(lng).words), keyboard, end
                        % keep numwords in order
                        nw=nw-1;
                        if ~isempty(nb) & any(nb==lng)
                            nwb = nwb-1;
                        end
                    end
                end
            end
        end

        % got to end of current interval
        %update n, nn, nnt, cnt
        cnt=nnt;
        % get rid of branching language from list of langs
        delni=find([s(n).time]==nnt);
        delsi=n(delni);
        n(delni)=[];
        % add in branching language's children
        newn=[s(delsi).child];
        n=[n,newn];

        if isempty(n)
            % no more langs in list - done
            break;
        else
            % update distance matrix
            if delsi > 2 & length([s(delsi).child]) > 2
                %XXX TODO change code below so that ancestral nodes introduced at the same time can borrow from each other
                warning('Polymorph may not be working when 2 ANST nodes are introduced at exactly the same time')
            end
            for si = delsi
                % copy row
                keeprow = dist(si,:);

                % set row and column to be infinite
                dist(si,:) = inf;
                dist(:,si) = inf;
                ci = s(si).child;
                for c = ci
                    dist(c,:) = keeprow;
                    dist(:,c) = keeprow';
                end
                if ~isempty(ci)
                    dist(ci(1),ci(2))= 0;
                    dist(ci(2),ci(1))= 0;
                end
            end
            nb = near(dist,maxdist);
            % allow cladegenic change
            for d=delsi
                for c=s(d).child
                    if p < 1
                        % have cladagenic loss
                        if min_wpm==0
                            [seq(c).words,seq(c).L]=BranchVocab(seq(d).words,p); %#ok<AGROW>
                        elseif min_wpm>=1
                            [seq(c).words,seq(c).L]=BranchVocabPolyMorph(seq(d).words,p,min_wpm);
                        else
                            error('Weirdness in PolyMorph(): min_wpm should be a whole number (for the moment, 0 or 1)');
                        end
                    else
                        % no cladagenic loss - just copy parent language
                        seq(c).words = seq(d).words;
                        seq(c).L = seq(d).L;
                         if seq(c).L~=numel(seq(c).words), keyboard, end
                    end
                end
            end
            % get next time
            nnt=max([s(n).time]);
        end
    end

    for k=1:(2*NS-1)
        if L_total==0, rseq{k}=[]; mseq{k}=[]; end
        if s(k).type>=ANST
            rseq{k}=[rseq{k},L_total+seq(k).words];
        elseif s(k).type==LEAF
            %if ~MISDAT || rand<s(k).xi % No missing data
                rseq{k}=[rseq{k},L_total+seq(k).words];
            %else % all data for this meaning category at this leaf go missing
            %    mseq{k}=[mseq{k},L_total+1:L_total+L];
            %end
        else
            disp('Warning: unknown node type in PolyMorph.m')
        end
        seq(k).words=[];
    end
    L_total=L_total+L;
    L_vals=[L_vals,L];
end
L=L_total;
%langs=num2cell(strjust(num2str([leaves]'),'left'),2);
langs = {s(leaves).Name}';
%keyboard;
cogs=num2cell(strjust(num2str([1:L_total]'),'left'),2);

function x = near(dist,maxdist)
x = find(sum(dist<maxdist) > 0);
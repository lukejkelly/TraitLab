function [rseq,langs,cogs,NS,L,ListOfCats]=SimBoInBlocks(br,mu,p,lambda,rho,kappa,nu,s,xi,allowempty,gomissing)
%[rseq,langs,cogs,NS,L,ListOfCats]=SimBo(br,mu,p,lambda,rho,kappa,nu,s)
%Simulate borrowing

global LEAF ROOT MIX

nmc=100; %number of meaning categories

leaves=find([s.type]==LEAF);
NS=length(leaves);
Root=find([s.type]==ROOT);

L=MyPoisson(lambda/mu);
seq(Root).words=1:L;
seq(Root).L=L;

cnt=s(Root).time;
n=[]; %the list of nodes we are currently looking at
for c=s(Root).child
   [seq(c).words,seq(c).L]=BranchVocab(seq(Root).words,p); %for p=1, this simply passes down all the vocabulary
   seq(c).cat=[];
end   
newn=[s(Root).child];
n=[n,newn];
nnt=max([s(n).time]);
ListOfCats=zeros(1,2*NS);
seq1=seq;
for k=1:nmc
    seqT=seq;
    while 1

       nw=sum([seqT(n).L]);
       numlang=length(n);

       y=0;
       T=cnt-nnt;
       while 1
          rate=nw*(br+mu)+numlang*(lambda+rho);
          y=y-log(rand)/rate; %time until next event
          if y>T, break; end; %go to next branching event
          r=rand;
          if r<nw*br/rate %borrowing event
             lngout=n(disample([seqT(n).L]./nw)); %choose 'out' language, weighted by number of words
             wd=ceil(rand*seqT(lngout).L);
             lngin=n(ceil(rand*numlang));
             nw=nw-seqT(lngin).L;
             seqT(lngin).words=union(seqT(lngin).words,seqT(lngout).words(wd));
             seqT(lngin).L=length(seqT(lngin).words);
             nw=nw+seqT(lngin).L;
          elseif r<nw*(br+mu)/rate %death event
             lng=n(disample([seqT(n).L]./nw));
             if allowempty || seqT(lng).L>0
                 wd=ceil(rand*seqT(lng).L);
                 seqT(lng).words(wd)=[];
                 seqT(lng).L=seqT(lng).L-1;
                 nw=nw-1;
             end
          else r<(nw*(br+mu) +numlang*lambda)/rate %birth event
             lng=n(ceil(rand*numlang));
             L=L+1;
             seqT(lng).words=[seqT(lng).words,L];
             seqT(lng).L=seqT(lng).L+1;
             nw=nw+1;
          else %catastrophic event
             lng=n(ceil(rand*numlang));
             seqT(lng).cat=[seqT(lng).cat, (cnt-y)];
             ListOfCats(lng)=ListOfCats(lng)+1;
             nw=nw-seqT(lng).L;
             seqT(lng).words=seqT(lng).words(rand(1,seqT(lng).L)>kappa);
             new=MyPoisson(nu);
             seqT(lng).words=[seqT(lng).words, (L+1):(L+new)];
             seqT(lng).L=length(seqT(lng).words);
             L=L+new;
             nw=nw+seqT(lng).L;
          end
       end

       %update n, nn, nnt, cnt
       cnt=nnt;
       delni=find([s(n).time]==nnt);
       delsi=n(delni);
       n(delni)=[];   
       newn=[s(delsi).child];
       n=[n,newn];
       if isempty(n)
          break;
       else
          for d=delsi
             for c=s(d).child
                [seqT(c).words,seqT(c).L]=BranchVocab(seqT(d).words,p);
                seqT(c).cat=[];
             end   
          end
          nnt=max([s(n).time]);
       end

    end
    
    for i=leaves
        if gomissing && rand>xi(i), seqT(i).words(:)=MIX; end
        seq(i).words=[seq(i).words seqT(i).words];
    end
end

[rseq{1:(2*NS-1)}]=deal(seq.words);
% for i=1:(2*NS-1)
%     TimeOfCats(i).times=seq(i).cat;
% end

%langs=num2cell(strjust(num2str([leaves]'),'left'),2);
langs = {s(leaves).Name}';
cogs=num2cell(strjust(num2str([1:L]'),'left'),2);

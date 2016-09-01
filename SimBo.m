function [rseq,langs,cogs,NS,L,ListOfCats,s]=SimBo(br,mu,p,lambda,rho,kappa,nu,s)
%[rseq,langs,cogs,NS,L,ListOfCats,s]=SimBo(br,mu,p,lambda,rho,kappa,nu,s)
%Simulate borrowing and catastrophes

global LEAF ROOT CHANGENU
if (exist('CHANGENU') & CHANGENU), nu=nu*CHANGENU; end

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


while 1
   
   nw=sum([seq(n).L]);
   numlang=length(n);
   
   y=0;
   T=cnt-nnt;
   while 1
      rate=nw*(br+mu)+numlang*(lambda+rho);
      y=y-log(rand)/rate; %time until next event
      if y>T, break; end; %go to next branching event
      r=rand;
      if r<nw*br/rate %borrowing event
          done=0;
          count=0;
          while ~done
             lngout=n(disample([seq(n).L]./nw)); %choose 'out' language, weighted by number of words
             wd=ceil(rand*seq(lngout).L);
             lngin=n(ceil(rand*numlang));
             if ~any(find(seq(lngin).words==seq(lngout).words(wd)))
                 nw=nw-seq(lngin).L;
                 seq(lngin).words=union(seq(lngin).words,seq(lngout).words(wd));
                 seq(lngin).L=length(seq(lngin).words);
                 nw=nw+seq(lngin).L;
                 done=1;
             else
                 count=count+1;
                 if count>50, disp('Could not find anything to borrow'); done=1;end
             end   
          end
             

      elseif r<nw*(br+mu)/rate %death event
         lng=n(disample([seq(n).L]./nw));
         wd=ceil(rand*seq(lng).L);
         seq(lng).words(wd)=[];
         seq(lng).L=seq(lng).L-1;
         nw=nw-1;
      elseif r<(nw*(br+mu) +numlang*lambda)/rate %birth event
         lng=n(ceil(rand*numlang));
         L=L+1;
         seq(lng).words=[seq(lng).words,L];
         seq(lng).L=seq(lng).L+1;
         nw=nw+1;
      else %catastrophic event
         lng=n(ceil(rand*numlang));
         seq(lng).cat=[seq(lng).cat, (cnt-y)];
         ListOfCats(lng)=ListOfCats(lng)+1;
         nw=nw-seq(lng).L;
         seq(lng).words=seq(lng).words(rand(1,seq(lng).L)>kappa);
         new=MyPoisson(nu);
         seq(lng).words=[seq(lng).words, (L+1):(L+new)];
         seq(lng).L=length(seq(lng).words);
         L=L+new;
         nw=nw+seq(lng).L;
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
            [seq(c).words,seq(c).L]=BranchVocab(seq(d).words,p);
            seq(c).cat=[];
         end   
      end
      nnt=max([s(n).time]);
   end
   
end

[rseq{1:(2*NS-1)}]=deal(seq.words);
for i=1:(2*NS-1)
    TimeOfCats(i).times=seq(i).cat;
end

%langs=num2cell(strjust(num2str([leaves]'),'left'),2);
langs = {s(leaves).Name}';
cogs=num2cell(strjust(num2str([1:L]'),'left'),2);

for i=1:2*NS-1
    s(i).cat=seq(i).cat;
end
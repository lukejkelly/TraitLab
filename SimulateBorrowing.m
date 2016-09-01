function [rseq,langs,cogs,NS,L]=SimulateBorrowing(br,mu,p,lambda,s)

%Written by Geoff Nicholls

GlobalSwitches
LEAF = 1;
ANST = 2;
ROOT = 3;
ADAM = 4;
leaves=find([s.type]==LEAF);
NS=length(leaves);
Root=find([s.type]==ROOT);
L=MyPoisson(lambda/mu);
seq(Root).words = 1:L;
seq(Root).L=L;
cnt=s(Root).time;
n=[];
for c=s(Root).child
   [seq(c).words,seq(c).L]=BranchVocab(seq(Root).words,p);
end  
newn=[s(Root).child];
n=[n,newn];
nnt=max([s(n).time]);
while 1
  
   nw=sum([seq(n).L]);
   numlang=length(n);
  
   y=0;
   T=cnt-nnt;
   while 1
      rate=nw*(br+mu)+numlang*lambda;
      y=y-log(rand)/rate;
      if y>T, break; end;
      r=rand;
      if r<nw*br/rate
         lngout=n(disample([seq(n).L]./nw));
         wd=ceil(rand*seq(lngout).L);
         %lngin=n(disample([seq(n).L]./nw));
         lngin=n(ceil(rand*numlang));
         nw=nw-seq(lngin).L;
         seq(lngin).words=union(seq(lngin).words,seq(lngout).words(wd));
         seq(lngin).L=length(seq(lngin).words);
         nw=nw+seq(lngin).L;
      elseif r<nw*(br+mu)/rate
         lng=n(disample([seq(n).L]./nw));
         wd=ceil(rand*seq(lng).L);
         seq(lng).words(wd)=[];
         seq(lng).L=seq(lng).L-1;
         nw=nw-1;
      else
         lng=n(ceil(rand*numlang));
         L=L+1;
         seq(lng).words=[seq(lng).words,L];
         seq(lng).L=seq(lng).L+1;
         nw=nw+1;
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
         end  
      end
      nnt=max([s(n).time]);
   end
  
end
[rseq{1:(2*NS-1)}]=deal(seq.words);
%langs=num2cell(strjust(num2str([leaves]'),'left'),2);
langs = {s(leaves).Name}';
cogs=num2cell(strjust(num2str([1:L]'),'left'),2);
function [x,y]=draw(s,f,TEXT,figlabel,cognate,cognomen,clade)

% [x,y]=draw(s,f,TEXT,figlabel,cognate,cognomen)


global LEAF ANST ROOT FONTSIZEINFIG COGS CONC MCMCCAT

% Luke 04/05/15 This statement does not have its intended purpose as after
% defining a global variable, it exists with the state it was defined to have
% elsewhere or is empty if it is being initialised for the first time. Either
% way, it always exists.
% if ~exist('MCMCCAT'), MCMCCAT=0; end
% Instead we use
if isempty(MCMCCAT), MCMCCAT = 0; end

if nargin == 1
   f=1;
   TEXT=0;
   figlabel='';
end

COGSUBTRE = 0;
if nargin>=5
   if ~isempty(cognate);
       COGSUBTRE=1;
   end
end

TL=size(s,2);
NS=TL/2;
N=TL-1;

LEAVES=[];
NODES=[];
Root=[];

for k=1:N
   if s(k).type==ROOT
       Root=k;
       NODES=[NODES,k];
   elseif s(k).type==LEAF
       LEAVES=[LEAVES,k];
   elseif s(k).type==ANST
       NODES=[NODES,k];
   end
end

nu=progeny(s,Root,LEAF);
n(nu(1,:))=nu(2,:);

v=split(s,n,0,1,Root);
x(v(1,:))=v(2,:);

figure(f);
clf(f);
currax=axes('Parent',f);

set(f,{'Name','NumberTitle'},{figlabel,'off'});
hold on;

for n=1:N
   y(n)=s(n).time;
   z1(n)=s(n).timedata(1);
   %z2(n)=s(n).timedata(2);
end

my=max(y);
y=my-y;
Tmax=my;
%z1=my-z1;
%z2=my-z2;

for k=NODES
   kc=s(k).child;
   for c=1:2
       clr='b';
       if COGSUBTRE
           if ismember(cognate,s(kc(c)).CovI)
               clr='r';
           elseif ismember(cognate,[s(kc(c)).ActI{:}])
               clr='k';
           end
       end
%       plot([z1(k),z2(k)],[x(k),x(k)],'g','Parent', currax,'LineWidth',2);
       plot([y(k),y(k)],[x(k),x(kc(c))],clr,'Parent', currax,'LineWidth',2);
       plot([y(k),y(kc(c))],[x(kc(c)),x(kc(c))],clr,'Parent', currax,'LineWidth',2);
   end
end

plot(y(Root),x(Root),'r*','Parent',currax);


%Plot catastrophes
for k=NODES
    if MCMCCAT && ~isempty(s(k).cat)
        for t=s(k).cat
            plot(Tmax-t,x(k),'ro','MarkerFaceColor','r','Parent',currax);
            
        end
    end
end

for k=LEAVES
    if MCMCCAT && ~isempty(s(k).cat)
        for t=s(k).cat
            plot(Tmax-t,x(k),'ro','MarkerFaceColor','r','Parent',currax);
        end
    end
end



yscale=(max(y)-min(y))/200;
xscale=(max(x)-min(x))/200;

if TEXT==LEAF
   for m=LEAVES
       text(y(m)+3*yscale,x(m),strrep(s(m).Name,'_','\_'),'Rotation',0,'FontSize',FONTSIZEINFIG,'Parent',currax);
   end
elseif TEXT==ANST
   for m=LEAVES
       text(y(m)+yscale,x(m),[num2str(m),'-',strrep(s(m).Name,'_','\_')],'Rotation',0,'FontSize',FONTSIZEINFIG,'Parent',currax);
   end
   for k=NODES
       text(y(k),x(k)-xscale,num2str(k),'Parent',currax);
   end
elseif TEXT==COGS
   for m=LEAVES
       text(y(m)+yscale,x(m),[num2str(m),'-',strrep(s(m).Name,'_','\_')],'Rotation',0,'FontSize',FONTSIZEINFIG,'Parent',currax);
   end
   for k=NODES
       text(y(k),x(k)-xscale,num2str(length(s(k).CovI)),'Parent',currax);
   end
 elseif TEXT==CONC
   for m=LEAVES
       text(y(m)+yscale,x(m),strrep(s(m).Name,'_','\_'),'Rotation',0,'FontSize',FONTSIZEINFIG,'Parent',currax);
   end
   for k=NODES
       text(y(k)+yscale,x(k),s(k).Name,'Parent',currax,'FontSize',FONTSIZEINFIG);
   end
end

%axis tight;
%axis fill;

xy=[min(y)-20*yscale,max(y)+20*yscale,min(x)-5*xscale,max(x)+5*xscale];
axis(currax,xy);
%set(currax, 'XTick', 0:500:max(get(currax, 'XTick'))) %%%%% REMOVE 17/09/2015
xy=axis(currax);
V=strjust(num2str(get(currax,'XTick')'),'left');
axis(currax, 'off');

%if ~(TEXT==CONC)
sV=size(V,1);
nV=max(y)-str2num(V);
yV=repmat(xy(3),sV,1);
ytV=repmat(xy(4),sV,1);
text(nV+max(y)*0.01,yV-(xy(4)-xy(3))/50,V,'Parent',currax);
plot([nV,nV]',[yV,ytV]','k:','Parent',currax);
%end

set(currax,'Position',[0 0.05 0.85 0.925]);

if COGSUBTRE
   if ~isempty(findstr(cognomen,'_'))
       text(nV(1),xy(4),['Cognate ',strrep(cognomen,'_','{_'),'}'],'Parent',currax);
   else
       text(nV(1),xy(4),['Cognate ',cognomen],'Parent',currax);
   end
end

if nargin==7,
   n=size(clade,2);
   for m=1:n,
       if size(clade{m}.rootrange,2)==2 & ~any(isinf(clade{m}.rootrange))
           U=GetLeaves(s,clade{m}.language);
           if ~isempty(U)
               i=mrca(U,s,Root);
               plot(y(i),x(i),'ks','MarkerFaceColor','k','MarkerSize',4);
               plot(my-clade{m}.rootrange,[x(i),x(i)]-3*xscale,'k-','LineWidth',2);
               plot([my-clade{m}.rootrange(1)].*[1,1],[x(i)-2*xscale,x(i)-4*xscale],'k-','LineWidth',2);
               plot([my-clade{m}.rootrange(2)].*[1,1],[x(i)-2*xscale,x(i)-4*xscale],'k-','LineWidth',2);
           end
       end
       if size(clade{m}.adamrange,2)==2 & ~any(isinf(clade{m}.adamrange))
           U=GetLeaves(s,clade{m}.language);
           if ~isempty(U)
               j=mrca(U,s,Root);
               i=s(j).parent;
               plot(y(i),x(i),'ms','MarkerFaceColor','m','MarkerSize',4);
               plot(my-clade{m}.adamrange,[x(i),x(i)]-3*xscale,'m-','LineWidth',2);
               plot([my-clade{m}.adamrange(1)].*[1,1],[x(i)-2*xscale,x(i)-4*xscale],'m-','LineWidth',2);
               plot([my-clade{m}.adamrange(2)].*[1,1],[x(i)-2*xscale,x(i)-4*xscale],'m-','LineWidth',2);
           end
       end
   end
end

if nargout==0, clear x,y; end

drawnow;

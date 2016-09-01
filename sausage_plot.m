function sausage_plot(x,y,col,F,eo,names)

figure(F);
r=8;
hold on;

sX=size(x);
len=sX(1);
dim=sX(2);
if (len==1) & (dim>1)
    x=x';
    y=y';
    len=sX(2);
    dim=sX(1);
end

count=0;
for loop=1:dim
    count=count+1;
    yp=y(:,loop);
    mv=max(yp);
    inds=find(yp>(mv/100));
    yp=yp./(r*mv);
    hold on;
    stairs(count+eo+yp(inds),x(inds),col{loop});
    stairs(count+eo-yp(inds),x(inds),col{loop});
end
set(gca,'XTick',[],'Xcolor',[1 1 1]);
set(gcf,'Color',[1 1 1]);

if nargin==6,
    xy=axis;
    text((1:dim),(xy(3)-0.02*(xy(4)-xy(3)))*ones(1,dim),names,'Rotation',270,'FontSize',8)
end

hold off;

function sausage_hist(V,col,F,eo,names)

figure(F);
r=6;

sV=size(V);
len=sV(1);
dim=sV(2);

count=0;
for loop=1:dim
    count=count+1;
    yp_vars=V(:,loop);
    max_yp=max(yp_vars);
    min_yp=min(yp_vars);
    range=max_yp-min_yp;
    N=sqrt(len);
    offset=(range)/(2*N);
    [a,b]=hist(yp_vars,linspace(min_yp+offset,min_yp+(range-offset),N));
%     tot=len; ma=max(a);
%     while tot>(0.95*len)
%         [i,j]=min(a);
%         a(j)=ma+1;
%         tot=tot-i;
%     end
%     a(a>ma)=0;
    mv=r*max(a);
    a=a./mv;
    if mv==0, keyboard; end
    b=b-offset;
    hold on;
    i=1;
    while i<length(b)
        while i<length(b) && a(i)==0
            i=i+1;
        end
        j=i+1;
        while j<=length(b) && a(j)>0
            j=j+1;
        end
        j=j-1;
        if (i<length(b) || a(end)>0 )
            v1=count+[0,a(i:j),0]+eo;
            v2=count-[0,a(i:j),0]+eo;
            w=[v1,v2(end:-1:1)];
            x=[b(i),b(i:j),b(j)];
            y=[x,x(end:-1:1)];
            fill(w,y,col{count},'EdgeColor',[1 1 1]);
            %plot(count-[0,a(i:j),0]+eo,[b(i),b(i:j),b(j)],col);
        end
        i=j+1;
    end
end
set(gca,'XTick',[],'Xcolor',[1 1 1]);
set(gcf,'Color',[1 1 1]);

if nargin==5,
    xy=axis;
    text((1:dim),(xy(3)+0.05*(xy(4)-xy(3)))*ones(1,dim),names,'Rotation',270,'FontSize',6)
end
hold off;

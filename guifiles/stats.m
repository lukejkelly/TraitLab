function [fbar,stdf,tauf]=stats(v,mx,StatNames,fig)
% [fbar,stdf,tauf]=stats(v,mx,StatNames,fig)

[s,N]=size(v);

switch nargin
    case 1
        error('not enough inputs in stats.m')
    case 2
        StatNames = {};
        fig = figure;
    case 3
        fig = figure;
end

if ~isempty(fig)
    figure(fig)
    set(fig,'NumberTitle','off','Name','Autocorrelations');
end

tauf = zeros(1,s);
fbar = zeros(1,s);
stdf = zeros(1,s);

for f=1:s
    w=v(f,1:N);
    r=zeros(1,mx+1);
    N=length(w);
    for k=0:mx
        CC=correl1([w(1:(N-k))',w((1+k):N)']);
        r(k+1)=CC(1,2);
    end
    
    G=zeros(1,mx);
    for k=1:mx
        G(k)=r(k+1)+r(k);
    end
    
    tauf(f)=-r(1)+2*G(1); %r(1)=1
    avt=1;
    for M=1:(mx-1)
        if G(M+1)<G(M) && G(M+1)>0
            tauf(f)=tauf(f)+2*G(M+1);
            avt=avt+2*r(M+1)^2;
        else
            break;
        end
    end
    if tauf(f)<1, tauf(f)=1; M=1; end
    %if M~=1, disp('unexpected M>1 in stats');keyboard;pause; end, end %when r(2)<0
    %GKN 24/8/10 there is no contradiction if tau<1 and M=2
    
    avt=avt/N;
    ast=sqrt(avt);
    
    fbar(f)=mean(w);
    stdf(f)=sqrt(tauf(f)/N)*std(w);
    
    if ~isempty(fig)
        figure(fig);
        hdl = subplot(s,2,2*f-1);
        plot(0:mx,r,[0,mx],[0,0],'r',[0,mx],[2*ast,2*ast],'r--',[0,mx],[-2*ast,-2*ast],'r--','Parent',hdl);
        set(hdl,'XLim',[0 mx]);
        if ~isempty(StatNames)
            title([StatNames{f},' \tau_f=',num2str(tauf(f)),' M=',num2str(M),' N=',num2str(N)],'Parent',hdl);
        else
            title([' \tau_f=',num2str(tauf(f)),' M=',num2str(M),' N=',num2str(N)],'Parent',hdl);
        end
        ylabel('autocorrelation','Parent',hdl);
        hdl = subplot(s,2,2*f);
        plot(w,'Parent',hdl);
        set(hdl,'XLim',[0 N]);
        if nargin==3
            ylabel(StatNames{f},'Parent',hdl);
        end
    end
end
if ~isempty(fig)
    figure(fig);
    hdl = subplot(s,2,2);
    title('MCMC output trace','Parent',hdl);
    hdl = subplot(s,2,2*s-1);
    xlabel('lag','Parent',hdl);
    hdl = subplot(s,2,2*s);
    xlabel('MCMC updates','Parent',hdl);
end

function [r,n] = correl1(x)
%CORREL1 Compute correlation matrix with error checking.

[n,m] = size(x);
c = cov(x);
d = sqrt(diag(c)); % sqrt first to avoid under/overflow
dd = d*d'; dd(1:m+1:end) = diag(c); % remove roundoff on diag
if any(dd==0)
    r = ones(size(dd));
else
    r = c ./ dd;
end
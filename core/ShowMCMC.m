function ShowMCMC(model,state,output,true)

%GlobalSwitches;
global GRAPH JUSTS JUSTT LEAF

N=output.Nsamp;
s=state.tree;
LEAVES=state.leaves;
Root=state.root;
NS=state.NS;

SYNTHETIC=(~isempty(true.state));
if SYNTHETIC,
  muTR=true.state.mu;
  pTR=true.state.p;
  LLtr=true.state.loglkd;
  LPtr=true.state.logprior;
  RootTR=true.state.root;
  sTR=true.state.tree;    
  LEAVEStr=true.state.leaves;
  s2str=[];
  for j=LEAVES
    leafname = s(j).Name;
    leafname = strrep(leafname,' ','');
    for k=LEAVEStr
      leafnameTR = sTR(k).Name;
      leafnameTR = strrep(leafnameTR,' ','');
      if isequal(leafname,leafnameTR)
        s2str=[s2str,k];
  	    break;
  	  end
    end
  end  
end

if any(output.verbose==[GRAPH JUSTS])

    xval=1:N;   
    figure(output.statsfig);
    set(output.statsfig,{'Name','NumberTitle'},{'MCMC output statistics','off'});
    hdl(1)=subplot(2,2,1); set(hdl(1),{'Parent'},{output.statsfig}); plot(xval,output.stats(1,:),'Parent',hdl(1));title('Log Prior','Parent',hdl(1));
    hdl(2)=subplot(2,2,2); set(hdl(2),{'Parent'},{output.statsfig}); plot(xval(ceil(N/4):N),output.stats(2,ceil(N/4):N),'Parent',hdl(2));title('Log Likelihood','Parent',hdl(2));
    hdl(3)=subplot(2,2,3); set(hdl(3),{'Parent'},{output.statsfig}); plot(xval,(output.stats(3,:)-min([s.time])),'Parent',hdl(3));title('Root Time','Parent',hdl(3));
    hdl(4)=subplot(2,2,4); set(hdl(4),{'Parent'},{output.statsfig}); plot(xval,LossRate(output.stats(4,:)),'Parent',hdl(4));title('Proportion Lost Traits','Parent',hdl(4));
    
    if SYNTHETIC,
        ERootTR=mrca(s2str,sTR,RootTR);
        mrcaT=sTR(ERootTR).time-min([sTR(s2str).time]);       
        set(hdl,'NextPlot','add');
        ylim =[min(LossRate([output.stats(4,:),muTR]))*0.9,max(LossRate([output.stats(4,:), muTR]))*1.1];
        %if ylim(2)<ylim(1), ylim=[0,1]; end
        set(hdl(4),'YLim',ylim);
        plot([xval(1),xval(N)],LossRate([muTR,muTR]),'Parent',hdl(4));
        plot([xval(1),xval(N)],[LPtr,LPtr],'Parent',hdl(1));
        plot([xval(ceil(N/4)),xval(N)],[LLtr,LLtr],'Parent',hdl(2));
        plot([xval(1),xval(N)],[mrcaT,mrcaT],'Parent',hdl(3));
        set(hdl,'NextPlot','replace');        
     end;
    
    drawnow;

end

if state.NS==2 & N>10 & any(output.verbose==[GRAPH JUSTS])

   set(hdl(3),'NextPlot','add');
   if SYNTHETIC
      k=common(s2str(1),s2str(2),sTR,RootTR);
      MRCA=sTR(k).time-min([sTR(s2str).time]);
      plot([xval(1),xval(N)],[MRCA,MRCA],'--','Parent',hdl(3));
      %nMLE=NaiveMLE(s(1).ActI{3},s(2).ActI{3},muTR);
      %plot([xval(1),xval(N)],[nMLE,nMLE],'r--','Parent',hdl(3));
   end      
   pM=pairMAP(s(1).ActI{3},s(2).ActI{3},state.mu,state.p);
   plot([xval(1),xval(N)],[pM,pM],'r--','Parent',hdl(3));
   set(hdl(3),'NextPlot','replace');
   drawnow;
   
   figure(output.postfig);
   clf(output.postfig);
   set(output.postfig,{'Name','NumberTitle'},{'root age, posterior distribution','off'});
   [a,b]=hist(output.stats(3,ceil(N/10):N)-min([s.time]),30);
   fx=linspace(0.5*min(b),max(b)*2,150);
   f=npd(model.prior,s(Root).ActI,state.mu,state.p,fx);
   postax = axes('Parent',output.postfig);
   f=exp(f-max(f));
   %hdl = bar(exp(b),a);
   hdl = bar(b,a);
   set(hdl,'Parent',postax);
   f=N*(b(2)-b(1))*f./(sum(f)*(fx(2)-fx(1)));
   set(postax,'NextPlot','add');
   %plot(exp(fx),f,'Parent',postax);
   plot(fx,f,'Parent',postax);
   plot([pM,pM],[0,max(a)],'r-*','Parent',postax);
   set(postax,'NextPlot','replace');
   drawnow;      
   
elseif state.NS~=2 & any(output.verbose==[JUSTT GRAPH])
    for i=1:2*state.NS
        if state.cat(i), s(i).cat=[1/(state.cat(i)+1):1/(state.cat(i)+1):1-1/(state.cat(i)+1)].*(s(s(i).parent).time-s(i).time)+s(i).time; else s(i).cat=[]; end
    end

   [x,y]=draw(s,output.treefig,LEAF,'MCMC state',[],[],model.prior.clade); 

end

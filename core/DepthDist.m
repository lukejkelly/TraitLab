function [x,y,ndata,content]=DepthDist(content,state,BLENDINONES)

%clear; 
%load outMC11; 
%[x,y]=DepthDist(data.content,state,1);
%figure;plot(x,y,'.');xy=axis; hold on; plot([xy(1),xy(2)],[xy(1),xy(2)]); hold off;

GlobalSwitches;

disp(sprintf('Begin Messages from DepthDist() ------------------------------*'));
disp('generating synthetic data');
%write the control variables into structures used by fullsetup
fsu=pop('fsu');
fsu.RUNLENGTH         = 1     ;
fsu.SUBSAMPLE         = 1     ;
fsu.SEEDRAND          = OFF   ;
fsu.SEED              = 2     ;
fsu.DATASOURCE        = BUILD ;
fsu.DATAFILE          = ''    ;
fsu.DATASYN           = OFF   ;
fsu.SYNTHSTYLE        = OLDTRE    ;
fsu.SYNTHTREFILE      = ''   ;
fsu.SYNTHTRE          = state.tree; %data.true.state.tree   ;
fsu.TREEPRIOR         = YULE  ;
fsu.ROOTMAX           = 0     ;
fsu.MCMCINITTREESTYLE = EXPSTART    ;
fsu.MCMCINITTREEFILE  = ''    ;
fsu.MCMCINITTREE      = ''    ;
fsu.MCMCINITMU        = 0.18  ;
fsu.MCMCINITP         = 1     ;
fsu.MCMCINITTHETA     = 1e-3  ;
fsu.VERBOSE           = JUSTT ;
fsu.OUTFILE           = ''    ;
fsu.OUTPATH           = ''    ;
fsu.LOSTONES          = OFF    ;
fsu.LOST              = 0    ;
fsu.LOSSRATE          = LossRate(state.mu);  
fsu.PSURVIVE          = 1     ;
fsu.BORROW            = OFF    ;
fsu.BORROWFRAC        = 0    ;
fsu.LOCALBORROW       = OFF;
fsu.MAXDIST           = 0;
fsu.POLYMORPH         = OFF    ;
fsu.NMEANINGCLASS     = 1   ; 
fsu.MASKING           = OFF    ;
fsu.DATAMASK          = []    ;
fsu.ISCOLMASK         = OFF   ;
fsu.COLUMNMASK        = []   ;
fsu.NUMSEQ            = state.NS;  %data.true.NS    ;
fsu.VOCABSIZE         = state.lambda/state.mu;
fsu.THETA             = 0;
fsu.ISCLADE           = OFF   ;
fsu.CLADE             = []    ;
fsu.STRONGCLADES      = []    ;
fsu.GUITRUE           = []    ;
fsu.GUICONTENT        = []    ;
fsu.MISDAT            = MISDAT ;
fsu.BETA              = state.beta;

ndata = fullsetup(fsu);

disp(sprintf('\nMasking excess taxa from the data (taxa absent in tree\n)'));
[a,ai]=setdiff(content.language,ndata.content.language);

if BLENDINONES
    disp(sprintf('\nRemoving solo traits (singleton columns) from loaded data (for dist matrix comp)\n'));
    LT=1;
    CM=[];
    content=ObserveData(content,ai,CM,LT);
    disp(sprintf('\nBlending synthetic rare traits into loaded data (for dist matrix comp)\n'));
    oc=find(sum(ndata.content.array)==1);
    od=ndata.content.array(:,oc);
    [snl,ni]=sort(ndata.content.language);
    [sdl,di]=sort(content.language);
    if ~isequal(snl,sdl)
        keyboard;
    end
    [sdi,dii]=sort(di);
    pa=[content.array(di,:),od(ni,:)];
    barray=pa(dii,:);
    content.array=barray;
    [content.NS,content.L]=size(barray);
else
    LT=0;
    disp(sprintf('\nRemoving non-data traits (empty columns) from loaded data\n'));
    CM=[];
    content=ObserveData(content,ai,CM,LT);
    disp(sprintf('\nRemoving non-data traits (empty columns) from synthetic data\n'));
    ndata.content = ObserveData(ndata.content,[],CM,LT);
end

x=zeros(state.NS,state.NS); 
for k=2:state.NS, 
    for j=1:k, 
        x(j,k)=pairMAP(find(content.array(j,:)),find(content.array(k,:)),state.mu,1);
    end
end
x=triu(x);


y=zeros(state.NS,state.NS); 
[tl{1:state.NS}]=deal(state.tree(state.leaves).Name);
[stl,ti]=sort(tl); 
[sdl,di]=sort(content.language);
zerotime=min([state.tree(state.leaves).time]);
for k=2:state.NS
    dk=di(k);
    tk=state.leaves(ti(k));
    for j=1:k
        dj=di(j);
        tj=state.leaves(ti(j));
        m=mrca([tj,tk],state.tree,state.root);
        t=(2*state.tree(m).time-state.tree(tj).time-state.tree(tk).time)/2;
        y(dj,dk)=t; 
        y(dk,dj)=t;
    end
end
y=triu(y);
if BLENDINONES
    LT=1;
    CM=[];
    disp(sprintf('\nRemoving solo data-traits (singleton columns) from synthetic data generated in DepthDist()\n'));
    ndata.content = ObserveData(ndata.content,[],CM,LT);
    disp(sprintf('\nRemoving synthetic solo data-traits added to loaded data by DepthDist()\n'));
    content = ObserveData(content,[],CM,LT);
end
disp(sprintf('End Messages from DepthDist() ------------------------------*'));
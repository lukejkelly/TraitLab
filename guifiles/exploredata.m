function exploredata(state,content,comp,dispmat,disphist,handles)

GlobalSwitches;

disp(sprintf('Begin Messages from exploredata() ------------------------------*'));

pro = 0;
if pro ==1;
profile on
end

if nargin<=5
    handles.makeHTMLet=0;
end


%TODO - make a checkbox in analgui.m to fix BLENDINONES
%YUK - if there are more than one percent (should be more like fifty
%percent) singleton columns and no zero columns - indicating the data
%is reasonably clean - then trust singleton columns                    
%if (mean(sum(content.array)==1)>0.01) & (sum(sum(content.array)==0)==0) %TODO - XXX fix this up
%   BLENDINONES=OFF;
%    disp(sprintf('\nTrusting the solo cognates in the loaded data (retaining singleton columns for consistency check)\n'));
%    disp(sprintf('\nYou may want to kill these singleton columns and check results are unchanged\n'));
%else
    BLENDINONES=ON;
%    disp(sprintf('\nSolo cognates in the loaded data (singleton columns) are either missing or untrustworthy\n'));
%    disp(sprintf('\nSynthetic singleton columns will be blended into the loaded data for the consistency check\n'));
%end

disp(sprintf('\nCalculating estimates of root times of all possible 2 leaf subtrees\nBe a little patient\n'))
[HEt,y,syndata,cutcontent]=DepthDist(content,state,BLENDINONES); 

NS=state.NS;
if comp   
    numplots = {2,2};
else
    numplots = {1,2};
end

if comp
    % calculate pairwise estimates of root times for synthetic data
    Et=zeros(NS,NS);
    for i=1:(NS-1)
        ActIi=find(syndata.content.array(i,:)==IN);
        for j=(i+1):NS
            ActIj=find(syndata.content.array(j,:)==IN);
            Et(i,j)=pairMAP(ActIi,ActIj,state.mu,state.p);
        end
    end
end




opt = pop('output');

if disphist
    figure(opt.statsfig);
    if comp
        titletag = ' - data in file';
    else 
        titletag = '';
    end
    set(opt.statsfig,{'NumberTitle','Name'},{'off','Histograms of data'});
    % for real data
    % make histogram of number of languages each cognate is present in
    mw=max(sum(cutcontent.array==IN));
    subplot(numplots{:},1); hist(sum(cutcontent.array==IN),0:mw);
    xlim([-1,mw+1])
    title(['Taxa per trait' titletag])
    % make histogram of number of cognates in each language
    subplot(numplots{:},2); hist(sum(cutcontent.array==IN,2));
    xlim1 = xlim(gca);
    ax1 = gca;
    title(['Traits per taxon' titletag])
    if comp
        % for synthetic data
        titletag = ' - synthetic data';
        % make histogram of number of languages each cognate is present in
        subplot(numplots{:},3); hist(sum(syndata.content.array==IN),0:mw);
        xlim([-1,mw+1])
        title(['Taxa per trait' titletag])
        % for synthetic data
        % make histogram of number of cognates in each language
        subplot(numplots{:},4); hist(sum(syndata.content.array==IN,2));
        title(['Traits per taxon' titletag])
        xlim2 = xlim(gca);
        xlim([min([xlim1(1),xlim2(1)]),max([xlim1(2),xlim2(2)])])
        axes(ax1)
        xlim([min([xlim1(1),xlim2(1)]),max([xlim1(2),xlim2(2)])])
    end
    
    if handles.makeHTMLet
        writehtml([handles.HTMLpathet handles.HTMLnameet],'<h2>Histograms of data</h2>',opt.statsfig);
    end
end

if dispmat
    % plot distance matrix for real data
    figure(opt.postfig);
    clf(opt.postfig);
    set(opt.postfig,{'NumberTitle','Name'},{'off','Distance matrix'});
    ax7 = axes('Parent',opt.postfig);
    imagesc(HEt,'Parent',ax7);
    axis equal;
    axis off;
    set(ax7,'Position',[0 0.05 0.9 0.925]);
    text(NS.*ones(1,NS)+1,1:NS,strrep(cutcontent.language,'_','\_'));
    if handles.makeHTMLet
        writehtml([handles.HTMLpathet handles.HTMLnameet],'<h2>Distance matrix</h2>',opt.postfig);
    end
end

dispdistdepth=1;
if dispdistdepth
    % plot distance depth relation
    figure(opt.distdepth);
    clf(opt.distdepth);
    set(opt.distdepth,{'NumberTitle','Name'},{'off','Distance Depth Relation'});
    ax8 = axes('Parent',opt.distdepth);
    plot(HEt,y,'.','Parent',ax8);
    xlabel('Maximum a posteriori estimate of TMRCA for pair')
    ylabel('TMRCA for pair in given tree')
    xy=axis(ax8); 
    set(ax8,{'NextPlot'},{'add'});
    plot([xy(1),xy(2)],[xy(1),xy(2)],'Parent',ax8);
    set(ax8,{'NextPlot'},{'replace'});
    if handles.makeHTMLet
        writehtml([handles.HTMLpathet handles.HTMLnameet],'<h2>Distance Depth Relation</h2>',opt.distdepth);
    end
end
disp(sprintf('\nData analysis complete\n'))
disp(sprintf('End Messages from exploredata() ------------------------------*'));

if pro ==1
profile report explorerep
end

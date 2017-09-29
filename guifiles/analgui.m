function varargout = analgui(varargin)

% ANALGUI Application M-file for analgui.fig
% FIG = ANALGUI launch analgui GUI.
% ANALGUI('callback_name', ...) invoke the named callback.
GlobalSwitches;
GlobalValues;

if isstruct(varargin{1}) % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);
    handles.data = varargin{1};
    handles.output = varargin{2};

    if ~isfield(handles,'llkdax')
        % the handle for the axes has mysteriously disappeared
        % so replace it
       handles.llkdax = findobj(fig,'Type','axes');
   end

   hdl = [handles.currpossl handles.currposet handles.viewtreebutt handles.movbeginet handles.movendet handles.moviebutt handles.zoombutt];
   if handles.output.Nsamp >= 1
       % there is output loaded
       set(hdl,'Enable','on');
       plot(handles.output.stats(3,:),'Color','Black','Parent',handles.llkdax);
       ylim=get(handles.llkdax,'Ylim');
       if handles.output.Nsamp==1
           xlim = [0.9 handles.output.Nsamp];
           slidestep = [1 0.1];
       else
           xlim = [1 handles.output.Nsamp];
           slidestep = [1/(xlim(2)-xlim(1)) 0.1];
       end
       set(handles.llkdax,'Xlim',xlim);
       handles.line=line([1 1],ylim,'Color','Red');
       set(handles.currpossl,{'Min','Max','Value','SliderStep'},{xlim(1),xlim(2),1,slidestep});
       set(handles.currposet,'String','1');
       set([handles.movbeginet handles.movendet],{'String'},{num2str(xlim(1));num2str(xlim(2))});
       set(handles.laget,'String',num2str(ceil(handles.output.Nsamp/5)));
       set(handles.nsamptxt, 'String',num2str(handles.output.Nsamp));
       set(handles.titletxt,'String',sprintf('MRCA time trace for samples %1.0f to %1.0f',xlim));
       if ~isempty(handles.data.file)
           set(handles.datafiletxt,'String',handles.data.file);
           set(handles.datadirtxt,'String',handles.data.path);
       end

       if ~(isempty(handles.output.file))
           set([handles.outdirtxt handles.outfiletxt],{'String'},{strtrunc(handles.output.path,19);strtrunc(handles.output.file,19)})
       end

   else
       set(hdl,'Enable','off');
       disp('No output loaded')
   end

   if ~isfield(handles,'makeHTMLet'), handles.makeHTMLet=0; end

    guidata(fig, handles);

  % Decrease font sizes if we're on a windows computer.
  decreaseFontSizesIfReq(handles)

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch thiserr
		disp(thiserr);
	end

end


function makeHTMLcb_Callback(h, eventdata, handles)
if get(h,'Value')==1
    handles.makeHTMLet=1;
else
    handles.makeHTMLet=0;
end

if handles.makeHTMLet && (~isfield(handles,'HTMLnameet') || isempty(handles.HTMLnameet))
    % use default name
    %set([handles.HTMLpathet
    %handles.HTMLnameet],{'String'},{'';datestr(now,'yyyymmddHHMMSS')});
    handles.HTMLnameet=[datestr(now,'yyyymmddHHMMSS') '.html'];
    handles.HTMLpathet=pwd;
    set(handles.HTMLfiletxt,'String',strtrunc(handles.HTMLnameet,55));
    set(handles.HTMLdirtxt,'String',strtrunc(handles.HTMLpathet,55));
end

guidata(h,handles);



% --- Executes on button press in HTMLnamebutt.
function HTMLnamebutt_Callback(h, eventdata, handles)


[HTMLnameet HTMLpathet] = uiputfile({'*.html';'*.*'},'Choose HTML file name');

if HTMLnameet==0 % User pressed cancel
    HTMLnameet='';
    HTMLpathet='';
end

%[msgstr, msgid] = lastwarn;

if handles.makeHTMLet && (~exist('HTMLnameet','var') || isempty(HTMLnameet)) && (~isfield(handles,'HTMLnameet') || isempty(handles.HTMLnameet))
    % use default name
    HTMLnameet=datestr(now,'yyyymmddHHMMSS');
    HTMLpathet=pwd;
end

if handles.makeHTMLet && strcmpi(HTMLnameet(end-4:end),'.html')
    % remove .html extension (to avoid it appearing twice)
    HTMLnameet=HTMLnameet(1:end-5);
end

handles.HTMLpathet=HTMLpathet;
handles.HTMLnameet=HTMLnameet;
if ~isempty(handles.HTMLnameet)
	disp(['HTML file ', [HTMLpathet, HTMLnameet], '.html selected'])
else
	disp('No HTML file selected')
end

set(handles.HTMLfiletxt,'String',strtrunc( handles.HTMLnameet,55));
set(handles.HTMLdirtxt,'String',strtrunc( handles.HTMLpathet,55));

if handles.makeHTMLet && isfield(handles,'data') && isfield(handles.data,'file') && ~isempty(handles.data.file)
    writehtml([handles.HTMLpathet handles.HTMLnameet],['Using data file ' handles.data.path handles.data.file '. <br/>']);
end

if handles.makeHTMLet && isfield(handles,'outfile') && ~isempty(handles.outfile)
    writehtml([handles.HTMLpathet handles.HTMLnameet],['Using output file ' handles.output.path handles.output.file '. <br/>']);
end

guidata(h,handles);






% --------------------------------------------------------------------
function  currposssl_Callback(h, eventdata, handles, varargin)

pos = round(get(h,'Value'));
set(h,'Value',pos);
set(handles.currposet,'String',sprintf('%1.0f',pos));
set(handles.line,'XData',[pos pos]);

% --------------------------------------------------------------------
function  analbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

vals = get([handles.distmatcb, handles.histscb, handles.compcb],'Value');
[dispmat,disphist,comp] = deal(vals{:});
% check that there is some data loaded
if isempty(handles.data.array)
   disp('Can''t analyse data as there is none loaded')
   return
end

if sum([dispmat,disphist,comp]) == 0
    disp('Select at least one option to do analysis')
else
    content = pop('content');
    content.array = handles.data.array;
    content.language = handles.data.language;
    [content.NS,content.L] = size(content.array);
    %TODO below would set up state properly but need model.prior

    %junk_model=pop('model');
    pos = str2double(get(handles.currposet,'string'));
    %state=makestate(junk_model.prior,handles.output.stats(4,pos),handles.output.stats(5,pos),content,rnextree(handles.output.trees{pos}));

    % check that there is some output loaded
    if pos == 0
        disp('Can''t analyse data as there is no output loaded')
        return
    end

    %step=ceil(handles.output.Nsamp/500); % DW removed unused variable
    state = pop('state');
    state.tree=rnextree(handles.output.trees{pos});
    state.mu = mean(handles.output.stats(4,:));
    state.p=mean(handles.output.stats(5,pos));
    state.lambda=mean(handles.output.stats(7,pos));
    state.root = find([state.tree.type]==ROOT);
    state.leaves = find([state.tree.type]==LEAF);
    state.nodes = [state.root, state.leaves, find([state.tree.type]==ANST)];
    state.NS = length(state.leaves);
    % check that have
%     state = pop('state');pos = str2double(get(handles.currposet,'string'));
%     state.tree = rnextree(handles.output.trees{pos});
%     state.mu = handles.output.stats(4,pos);
%     state.p=handles.output.stats(5,pos);
%     state.lambda=handles.output.stats(7,pos);
%     state.root = find([state.tree.type]==ROOT);
%     state.leaves = find([state.tree.type]==LEAF);
%     state.nodes = [state.root, state.leaves, find([state.tree.type]==ANST)];
%     state.NS = length(state.leaves);
    exploredata(state,content,comp,dispmat,disphist,handles); %GKN 7/7 dropped useless filename and HEt args
    guidata(h,handles);
end



% --------------------------------------------------------------------
function  matrixbutt_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function  viewtreebutt_Callback(h, eventdata, handles, varargin)
GlobalSwitches;

pos = str2double(get(handles.currposet,'string'));
s = rnextree(handles.output.trees{pos});
scat= rnextree(handles.output.cattrees{pos});
s=CatTreeToList(scat,s);

if ~isempty(handles.data.clade)
    draw(s,handles.output.treefig,0,['Sample number ' num2str(pos)],[],[],handles.data.clade);
else
    draw(s,handles.output.treefig,0,['Sample number ' num2str(pos)]);
end
if handles.makeHTMLet
    writehtml([handles.HTMLpathet handles.HTMLnameet],['<h2>Sample number ' num2str(pos) ' </h2>'],handles.output.treefig);
end

% --------------------------------------------------------------------
function  currposet_Callback(h, eventdata, handles, varargin)

% check that the input is valid
entry = str2double(get(h,'string'));
lims = get(handles.currpossl,{'Min','Max'});
lims = [lims{:}];
if isnan(entry) || floor(abs(entry)) ~= entry || entry<lims(1) || entry > lims(2)
    errordlg(sprintf('Position must be a positive integer between %1.0f and %1.0f',lims(1), lims(2)),'Invalid Input','modal')
    % return to same value as slider
    set(h,'String',num2str(get(handles.currpossl,'Value')));
else
    % good input
    % set slider pos
    set(handles.currpossl,'Value',entry);
    % set line pos
    set(handles.line,'XData',[entry entry]);
end

% --------------------------------------------------------------------
function  moviebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

% minimum gap between frames in seconds
gap = 1/3;
% find start and end of movie
lims = str2double(get([handles.movbeginet handles.movendet],'String'));
if abs(diff(lims))>=1
    % there is more than one sample
    if lims(1)>lims(2)
        inc = -1;
    else
        inc = 1;
    end

    xlim = sort(lims);
    currxlim = get(handles.llkdax,'Xlim');
    if xlim(1)<currxlim(1) || xlim(2) > currxlim(2)
        % want to show movie outside current llkd plot - zoom in on movie region
        plot(xlim(1):xlim(2),handles.output.stats(2,xlim(1):xlim(2)),'Color','Black','Parent',handles.llkdax);
        ylim=get(handles.llkdax,'Ylim');
        set(handles.llkdax,'Xlim',xlim);
        handles.line=line([1 1],ylim,'Color','Red');
        slidestep = [1/(xlim(2)-xlim(1)) 0.1];
        set(handles.currpossl,{'Min','Max','Value','SliderStep'},{xlim(1),xlim(2),xlim(1),slidestep});
        set(handles.currposet,'String',xlim(1));
        set(handles.titletxt,'String',sprintf('Log Likelihood Trace for samples %1.0f to %1.0f',xlim));
        guidata(h,handles);
    end

    for pos = lims(1):inc:lims(2)
        % check for interuption
        if pos > lims(1) && ~ishandle(handles.output.treefig)
            disp('Movie interrupted')
            return
        end
        tic
        s = rnextree(handles.output.trees{pos});
        scat= rnextree(handles.output.cattrees{pos});
        s=CatTreeToList(scat,s);
        draw(s,handles.output.treefig,0,['Sample number ' num2str(pos)]);
        % keep currpos counters and line up with play
        % set slider pos
        set(handles.currpossl,'Value',pos);
        % set line pos
        set(handles.line,'XData',[pos pos]);
        % set counter text
        set(handles.currposet,'String',sprintf('%1.0f',pos));
        % make sure that enough time has passed between frames
        elapsed=toc;
        if elapsed<gap
            pause(gap-elapsed);
        end
    end
else
    disp('There must be at least two samples in the interval to show movie')
end

% --------------------------------------------------------------------
function  movbeginet_Callback(h, eventdata, handles, varargin)
% check that the input is valid
checknum(h,1,handles.output.Nsamp,1,1);

% --------------------------------------------------------------------
function  movendet_Callback(h, eventdata, handles, varargin)
% check that the input is valid
checknum(h,1,handles.output.Nsamp,handles.output.Nsamp,1);

% --------------------------------------------------------------------
function  autocorrbutt_Callback(h, eventdata, handles, varargin)

maxlag = str2double(get(handles.laget,'string'));
first = str2double(get(handles.starthistet,'string'));
if maxlag > size(handles.output.stats,2) - first - 4
    fprintf('Value of Lag should be much smaller than Number of Samples being analysed (currently %g)\n',size(handles.output.stats,2) - first)
    return
end
vals = get([handles.priorcorrcb, handles.llkdcorrcb, handles.rootcorrcb, handles.mucorrcb handles.kappacorrcb handles.rhocorrcb handles.betacorrcb handles.betamucorrcb],'Value'); % LUKE 04/09/2016
vals = find([vals{:}]);
valsindex=[1 2 3 4 8 9 12 13]; % LJK

% Computing beta / mu, if necessary. LUKE 04/09/2016
if handles.betamucorrcb.Value;
   handles.output.stats(13, :) = handles.output.stats(12, :) ./ handles.output.stats(4, :);
end

if ~isempty(vals)
    names = {'PRIOR','LIKELIHOOD','ROOT TIME','MU','KAPPA','RHO', 'BETA', 'BETA / MU'};
    stats(handles.output.stats(valsindex(vals),first:end),maxlag,names(vals),handles.output.statsfig);
    if handles.makeHTMLet, writehtml([handles.HTMLpathet handles.HTMLfileet], '<h2>Autocorrelations</h2>',handles.output.statsfig); end
else
	fprintf('No statistics selected for which to calculate autocorrelations\n')
end

% --------------------------------------------------------------------
function  laget_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry) || floor(abs(entry)) ~= entry || entry > handles.output.Nsamp - 2
    errordlg(sprintf('Maximum lag must be a positive integer between 1 and %1.0f',handles.output.Nsamp - 2),'Invalid Input','modal')
    set(h,'String',num2str(ceil(handles.output.Nsamp/5)));
end

% --------------------------------------------------------------------
function  leafnamesrb_Callback(h, eventdata, handles, varargin)
set([handles.anstnumsrb handles.showcovrb],'Value',0);

% --------------------------------------------------------------------
function  anstnumsrb_Callback(h, eventdata, handles, varargin)

set([handles.leafnamesrb handles.showcovrb],'Value',0);

% --------------------------------------------------------------------
function  showcovrb_Callback(h, eventdata, handles, varargin)
set([handles.anstnumsrb handles.leafnamesrb],'Value',0);


% --------------------------------------------------------------------
function  showcogcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==1
    set([handles.showcoget],'Enable','on');
else
    set([handles.showcoget],'Enable','off');
end


% --------------------------------------------------------------------
function  showcoget_Callback(h, eventdata, handles, varargin)
% check that the input is valid
entry = str2double(get(h,'String'));
L = size(handles.data.array,2);
if isnan(entry) || floor(abs(entry)) ~= entry || entry > L
    errordlg(sprintf('Start must be a positive integer between 1 and %1.0f',L),'Invalid Input','modal')
    set(h,'String','1');
end

% --------------------------------------------------------------------
function  loadoutbutt_Callback(h, eventdata, handles, varargin)

[filename pathname] = uigetfile({'*.nex','Nexus'},'Load an output file');

if isequal(filename,0)||isequal(pathname,0)
    %no file selected - done
    disp('No file selected')
    return
else
    %user selected file - check that it is a .txt or a .nex file
    if strcmp(strtrunc(filename,4),'.nex') || strcmp(strtrunc(filename,4),'.txt')
        % good file - load output
        disp(['Output file ', [pathname, filename(1:end-4)], ' selected'])
        handles.output=readoutput([pathname, filename(1:end-4)]);
        set([handles.outdirtxt handles.outfiletxt],{'String'},{pathname;filename});
        handles.output.path=pathname;
        handles.output.file=filename(1:end-4);
        if handles.makeHTMLet, writehtml([handles.HTMLpathet handles.HTMLnameet],['Using output file ' handles.output.path handles.output.file '. <br/>']); end
        guidata(h,handles);
        if handles.output.Nsamp >= 1
            % there is output loaded
            hdl = [handles.currpossl handles.currposet handles.viewtreebutt handles.movbeginet handles.movendet handles.moviebutt handles.zoombutt];
            set(hdl,'Enable','on');
            plot(handles.output.stats(2,:),'Color','Black','Parent',handles.llkdax);
            ylim=get(handles.llkdax,'Ylim');
            if handles.output.Nsamp==1
                xlim = [0.9 handles.output.Nsamp];
                slidestep = [1 0.1];
            else
                xlim = [1 handles.output.Nsamp];
                slidestep = [1/(handles.output.Nsamp-1) 0.1];
            end
            set(handles.llkdax,'Xlim',xlim);
            handles.line=line([1 1],ylim,'Color','Red');
            set(handles.currpossl,{'Min','Max','Value','SliderStep'},{xlim(1),xlim(2),1,slidestep});
            set(handles.currposet,'String','1');
            set([handles.movbeginet handles.movendet],{'String'},{num2str(xlim(1));num2str(xlim(2))});
            set(handles.laget,'String',num2str(ceil(handles.output.Nsamp/5)));
            set(handles.nsamptxt, 'String',num2str(handles.output.Nsamp));
            set(handles.titletxt,'String',sprintf('Log Likelihood Trace for samples %1.0f to %1.0f',xlim));
            guidata(h,handles);
        else
            hdl = [handles.currpossl handles.currposet handles.viewtreebutt handles.movbeginet handles.movendet handles.moviebutt handles.zoombutt];
            set(hdl,'Enable','off');
            set(handles.titletxt,'String','LOG LIKELIHOOD PLOT');
            cla(handles.llkdax);
            set(handles.nsamptxt, 'String','0');
            disp('No output loaded')
        end

    else
        % bad file type
        disp('No file opened - file type must be either .nex or .txt')
    end
end

% --------------------------------------------------------------------
function  insptreebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
global za zb zc za1 za2 zb1 zb2 zq1 zq2
% merge data and tree
pos = str2double(get(handles.currposet,'string'));
s = rnextree(handles.output.trees{pos});
scat= rnextree(handles.output.cattrees{pos});
s=CatTreeToList(scat,s);
[s,ok]=mergetreedata(s,handles.data.array,handles.data.language);
if ok
    % merge worked, now need to fill work variables in tree
    state = pop('state');
    state.NS = length(s)/2;
    state.L = size(handles.data.array,2);
    state.tree =s;
    state.root = find([s.type] == ROOT);
    state.leaves = find([s.type] == LEAF);
    state.nodes = find([s.type]==ANST);
    [state.tree([state.leaves state.nodes]).mark]=deal(1);
    za=zeros(1,state.L);
    zb=zeros(1,state.L);
    zc=zeros(1,state.L);
    za1=zeros(1,state.L);
    zb1=zeros(1,state.L);
    zq1=zeros(1,state.L);
    za2=zeros(1,state.L);
    zb2=zeros(1,state.L);
    zq2=zeros(1,state.L);

   [state.tree]=ActiveI(state.tree,state.root);

   state.tree(state.root).CovI=[state.tree(state.root).ActI{:}];

   state.tree=CoversI(state.tree,state.root);

    % get drawing options
    vals = get([handles.leafnamesrb handles.anstnumsrb handles.showcovrb handles.showcogcb],'Value');

    if vals{1}
        vb = LEAF;
    elseif vals{2}
        vb = ANST;
    elseif vals{3}
        vb = COGS;
    else
        vb = -1;
    end
    showcog = 0;
    cogname = '';

    if vals{4}
        showcog = str2double(get(handles.showcoget,'String'));
        if length(handles.data.cognate)>=showcog
            cogname = handles.data.cognate{showcog};
        end
    end
    draw(state.tree,handles.output.treefig,vb,['Sample number ' num2str(pos)],showcog,cogname);


    if handles.makeHTMLet
        if vb == LEAF
            showing='Showing leaf names.';
        elseif vb==ANST
            showing='Showing all names/numbers.';
        elseif vb== COGS
            showing='Showing number of active traits.';
        else
            showing='';
        end
        if vals{4}, showing=[showing '<br/> Showing trait number' int2str(showcog) '.']; end
        writehtml([handles.HTMLpathet handles.HTMLnameet], ['<h2>Tree sample number ' num2str(pos) '</h2>' showing '<br/>'], handles.output.treefig);
    end
else
    disp('Unable to inspect tree due to incompatibility in loaded data and output')
end

% --------------------------------------------------------------------
function  priorcorrcb_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function  llkdcorrcb_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function  rootcorrcb_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function  mucorrcb_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function  betacorrcb_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function  betamucorrcb_Callback(h, eventdata, handles, varargin)
% --------------------------------------------------------------------
function  zoombutt_Callback(h, eventdata, handles, varargin)

xlim = sort(str2double(get([handles.movbeginet handles.movendet],'String')));
if xlim(1)<xlim(2)
    plot(xlim(1):xlim(2),handles.output.stats(2,xlim(1):xlim(2)),'Color','Black','Parent',handles.llkdax);
    ylim=get(handles.llkdax,'Ylim');
    set(handles.llkdax,'Xlim',xlim);
    handles.line=line([1 1],ylim,'Color','Red');
    slidestep = [1/(xlim(2)-xlim(1)) 0.1];
    set(handles.currpossl,{'Min','Max','Value','SliderStep'},{xlim(1),xlim(2),xlim(1),slidestep});
    set(handles.currposet,'String',xlim(1));
    set(handles.titletxt,'String',sprintf('Log Likelihood Trace for samples %1.0f to %1.0f',xlim));
    guidata(h,handles);
else
    disp('There must be at least two samples in the interval to zoom')
end



% --------------------------------------------------------------------
function  loaddatabutt_Callback(h, eventdata, handles, varargin)
% get file names from user
[filename,pathname] = uigetfile({'*.nex','Nexus'},'Choose a nexus file to load data from');
if isequal(filename,0)||isequal(pathname,0)
    %no file selected - done
    disp('No file selected')
else
    %user selected file - save info in relevant place
    disp(['Data file ', pathname, filename, ' selected'])
    handles.data.file = filename(1:end-4);
    handles.data.path = pathname;
    [s,content,true,clade] = nexus2stype([pathname filename]);
    handles.data.array = content.array;
    handles.data.language = content.language;
    handles.data.cognate = content.cognate;
    handles.data.clade=clade;
    set([handles.datadirtxt,handles.datafiletxt],{'String'},{pathname;filename});
    if handles.makeHTMLet, writehtml([handles.HTMLpathet handles.HTMLnameet],['Using data file ' handles.data.path handles.data.file '.nex. <br/>']); end
    guidata(h,handles);
end



% --------------------------------------------------------------------
function  distmatcb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function  histscb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')
    able(handles.compcb,[],handles.compcb,[]);
else
    able([],handles.compcb,[],handles.compcb);
end


% --------------------------------------------------------------------
function  compcb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function  histbutt_Callback(h, eventdata, handles, varargin)
bins = str2double(get(handles.binset,'string'));
starthist = str2double(get(handles.starthistet,'string'));
vals = get([handles.priorcorrcb, handles.llkdcorrcb, handles.rootcorrcb, handles.mucorrcb, handles.kappacorrcb, handles.rhocorrcb, handles.betacorrcb, handles.betamucorrcb],'Value'); % LUKE 04/09/2016
vals = find([vals{:}]);
valsindex=[1 2 3 4 8 9, 12, 13]; % LJK 12, 13

% Computing beta / mu, if necessary. LUKE 04/09/2016
if handles.betamucorrcb.Value;
   handles.output.stats(13, :) = handles.output.stats(12, :) ./ handles.output.stats(4, :);
end

if ~isempty(vals)
    names = {'PRIOR','LIKELIHOOD','ROOT TIME','MU','KAPPA','RHO', 'BETA', 'BETA / MU'};
    figure(handles.output.histfig);
    set(handles.output.histfig,'NumberTitle','off','Name','Histograms of output trace');
    for i = 1:length(vals)
        subplot(1,length(vals),i);
        hist(handles.output.stats(valsindex(vals(i)),starthist:end),bins);
        title(names(vals(i)));
        fprintf([names{vals(i)},' Mean %g, Std %g \n'],mean(handles.output.stats(valsindex(vals(i)),starthist:end)),std(handles.output.stats(valsindex(vals(i)),starthist:end)));
        [hpdl hpdh]=hpd(handles.output.stats(valsindex(vals(i)),starthist:end),5);
        fprintf('95%% HPD interval (assumes unimodal): [%g,%g] \n',hpdl,hpdh);
    end
    if handles.makeHTMLet
        writehtml([handles.HTMLpathet handles.HTMLnameet],'<h2>Histograms of output trace</h2>',handles.output.histfig);
    end
else
    fprintf('No statistics selected for which to make histograms\n')
end


% --------------------------------------------------------------------
function  binset_Callback(h, eventdata, handles, varargin)
checknum(handles.binset,2,100,20,1);


% --------------------------------------------------------------------
function  starthistet_Callback(h, eventdata, handles, varargin)
checknum(handles.binset,1,handles.output.Nsamp,1,1);




% --------------------------------------------------------------------
function  agehistbutt_Callback(h, eventdata, handles, varargin)
global ROOT LEAF

x = getintegertext(handles.agecompet,'Languages');
if ~isempty(x)
    % check that x is a vector of numbers of language tips
    pos = str2double(get(handles.currposet,'string'));
    pos_s = rnextree(handles.output.trees{pos});
    leaves = find([pos_s.type]==LEAF);
    if length(union(x,leaves))>length(leaves)
        disp('The vector of taxon numbers doesn''t match the numbers on tree leaves.')
        disp('Note that taxon numbers may change when the current tree changes.')
        disp('Use following list or Inspect Tree button to find current taxon numbers.')
        disp('Taxon number-Taxon name')
        for i = 1:length(leaves)
            fprintf('%1g - %1s \n',leaves(i),pos_s(leaves(i)).Name)
        end
    else
        langs={pos_s(x).Name};
        langlist=strrep(['(',sprintf('%s, ',langs{1:(end-1)}), sprintf('%s)',langs{end})],'_','\_');
        disp([sprintf('Computing MRCA distribution for '),langlist]);
        starthist = str2double(get(handles.starthistet,'string'));
        bins = str2double(get(handles.binset,'string'));
        maxlag = str2double(get(handles.laget,'string'));
        if maxlag > size(handles.output.stats,2) - starthist - 4
            fprintf('Value of Lag should be much smaller than Number of Samples being analysed (currently %g)\n',size(handles.output.stats,2) - first)
            return
        end


        if handles.makeHTMLet
            writehtml([handles.HTMLpathet handles.HTMLnameet],['<h2>MRCA distribution</h2><b2/>' int2str(bins) 'bins; ' 'lag: ' int2str(maxlag) '<br/>']);
        end

        for i = starthist:length(handles.output.trees)
            s = rnextree(handles.output.trees{i});
            sroot = find([s.type]==ROOT);
            y=[];
            for k=1:length(x),
                fi=find(strcmp(langs(k),{s.Name}));
                if length(fi)~=1, error('One of the taxa was not found or was not unique'); end
                y=[y,fi]; %#ok<AGROW>
            end
            mrcatime(i-starthist+1) = s(mrca(y,s,sroot)).time; %#ok<AGROW>
            cladeindicator(i-starthist+1) = IsClade(y,s,sroot); %#ok<AGROW>
        end
        disp(['Tree, TMRCA for ' langlist])
        %disp(sprintf('%g, %g \n',[starthist:length(handles.output.trees);mrcatime]))
        figure(handles.output.tmrcafig)
        set(handles.output.tmrcafig,{'Name','NumberTitle'},{'TMRCA distribution','off'});
        hist(mrcatime,bins);title(['MRCA-time distribution for clade ',langlist]);
        if handles.makeHTMLet, writehtml([handles.HTMLpathet handles.HTMLnameet], ['TMRCA for ' langlist '<br/>'],handles.output.tmrcafig); end

        stats(mrcatime,maxlag,{['MRCA: ',langlist]},handles.output.histfig);
        mmr=mean(mrcatime); smr=std(mrcatime);
        fprintf('mean MRCA-time: %g (std err %g) or [%g,%g] at 2*sigma \n',mmr,smr,mmr-2*smr,mmr+2*smr);
        [hpdl hpdh]=hpd(mrcatime,5);
        fprintf('95%% HPD interval (assumes unimodal): [%g,%g] \n',hpdl,hpdh);
        if handles.makeHTMLet
            writehtml([handles.HTMLpathet handles.HTMLnameet], sprintf('mean MRCA-time: %g (std err %g) or [%g,%g] at 2*sigma <br/>',mmr,smr,mmr-2*smr,mmr+2*smr));
            writehtml([handles.HTMLpathet handles.HTMLnameet], sprintf('95&#37; HPD interval (assumes unimodal): [%g,%g] <br/>',hpdl,hpdh),handles.output.histfig);
        end

        [fbar,stdf]=stats(cladeindicator,maxlag,{['IsClade(',langlist,')']},handles.output.histfig2);
        mci=mean(cladeindicator); %sci=std(cladeindicator);

        disp(['probability that ',langlist,sprintf(' is a clade: %g (std OF MEAN %g)',mci,stdf)]);
        if handles.makeHTMLet, writehtml([handles.HTMLpathet handles.HTMLnameet], ['probability that ',langlist,sprintf(' is a clade: %g (std OF MEAN %g)',mci,stdf)],handles.output.histfig2); end
    end
end
% --------------------------------------------------------------------
function  agecompet_Callback(h, eventdata, handles, varargin)
getintegertext(h,'taxa');


% --- Executes on button press in createconsensusbutt.
function createconsensusbutt_Callback(hObject, eventdata, handles)
    GlobalSwitches
    starthist = str2double(get(handles.starthistet,'string'));
    p=str2double(get(handles.consensuspet,'String'));
    step=str2double(get(handles.consensusstepet,'String'));
    includecats=get(handles.consensuscatcb,'Value');

    if includecats
         consensuscat(strcat(handles.output.path,handles.output.file),p,'',step,starthist);
%         [scon,conclade]=consensuscat(strcat(handles.output.path,handles.output.file),p,'',step,starthist);
    else
        temp=MCMCCAT; %#ok<NODEF>
        MCMCCAT=0; %#ok<NASGU>
        scon=consensus(strcat(handles.output.path,handles.output.file,'.nex'),p,'',step,starthist);
%        [scon,clade]=consensus(strcat(handles.output.path,handles.output.file,'.nex'),p,'',step,starthist);
        fig = pop('output');
        draw(scon,fig.constree,CONC,'Consensus tree');
        MCMCCAT=temp; %#ok<NASGU>
    end
    if handles.makeHTMLet
        writehtml([handles.HTMLpathet handles.HTMLnameet], sprintf('<h2>Consensus tree</h2> Displaying nodes with at least %g&#37; support.<br/>',100*p),fig.constree);
    end



function consensusstepet_Callback(hObject, eventdata, handles)


if (handles.output.Nsamp == 0)
	checknum(handles.consensusstepet,1,100,1,1);
else
	checknum(handles.consensusstepet,1,handles.output.Nsamp,1,1);
end


function consensuspet_Callback(hObject, eventdata, handles)
checknum(handles.consensuspet,.5,1,.5,0);



% --- Call this when opening figure. Posted by a guy called Wouter at
% https://stackoverflow.com/questions/19843040/matlab-gui-compatibility-between-mac-and-windows-display
function decreaseFontSizesIfReq(handles)
% make all fonts smaller on a Windows computer
if ispc()
  for afield = fieldnames(handles)'
    afield = afield{1}; %#ok<FXSET>
    try %#ok<TRYNC>
      set(handles.(afield), 'FontSize', 8.5); % decrease font size
    end
  end
end

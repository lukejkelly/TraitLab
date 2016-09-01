function varargout = treegui(varargin)
% TREEGUI Application M-file for treegui.fig

GlobalSwitches;
GlobalValues;

if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
    
    % Use system color scheme for figure:
    set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
    
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    % initialise other data stored in handles
    handles.data  = struct('array',{[]},...
        'language',{''},...
        'cognate',{''},...
        'file',{''},'path',{''});
    handles.tree = struct('tree',{TreeNode([],[],[],[],[],[])},...
        'file',{''},'path',{''});
    handles.output =pop('output');
    handles.output.file='/TreeMCMCout';
    handles.output.path=cd;
    guidata(fig, handles);
    
    set([handles.outdirtxt handles.outfiletxt],{'String'},{strtrunc(handles.output.path,19);handles.output.file});
    set([handles.syndataoutdirtxt handles.syndataoutfiletxt],{'String'},{[cd '/'];'synthdata.nex'});

    
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
    catch
        disp(lasterr);
    end
    
end

% --------------------------------------------------------------------
function varargout = pausebutt_Callback(h, eventdata, handles, varargin)

if strcmp(get(h,'String'),'Pause')
    set(h,'String','Resume');
    set(handles.statustxt,'String','Paused');
    set(handles.startbutt,'Enable','off');
    writeoutput(handles.output);
    waitfor(h,'String','Pause')
else
    % in Resume mode
    set(h,'String','Pause');   
    set(handles.statustxt,'String','Running');
end

% --------------------------------------------------------------------
function varargout = stopbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
GlobalValues;

set(handles.startbutt,{'Enable','UserData'},{'on',~STOPRUN});
set([handles.pausebutt,h],'Enable','off'); 
set(handles.pausebutt,'String','Pause');    
set(handles.statustxt,'String','Idle');

% --------------------------------------------------------------------
function varargout = fixmurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
set([handles.randommurb handles.specmurb],'Value',0);
set([handles.muvalfixet,handles.muvalet],{'Enable'},{'on';'off'});

% --------------------------------------------------------------------
function varargout = randommurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
set([handles.fixmurb handles.specmurb],'Value',0);
set([handles.muvalfixet,handles.muvalet],{'Enable'},{'off'});


% --------------------------------------------------------------------
function varargout = muvalet_Callback(h, eventdata, handles, varargin)
entry = str2double(get(h,'string'));
if isnan(entry) | 0 >= entry | entry >= 1  
    errordlg('Mu must be a number between 0 and 1','Invalid Input','modal')
    set(h,'String',num2str(0.18));
end

% --------------------------------------------------------------------
function varargout = muvalfixet_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry) | entry <= 0 | entry >= 1  
    errordlg('Mu must be a number between 0 and 1','Invalid Input','modal')
    set(h,'String',num2str(0.18));
end


% --------------------------------------------------------------------
function varargout = randtreerb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
set([handles.spectreerb,handles.truetreerb],'Value',0);
set(handles.initthetaet,'Enable','on');
% --------------------------------------------------------------------
function varargout = spectreerb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
set([handles.randtreerb,handles.truetreerb],{'Value'},{0});
set(handles.initthetaet,'Enable','off');


% --------------------------------------------------------------------
function varargout = runet_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry) | floor(abs(entry)) ~= entry | entry==0
    errordlg('Run length must be a positive integer','Invalid Input','modal')
    set(h,'String',num2str(1e5));
end


% --------------------------------------------------------------------
function varargout = sampleet_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry) | floor(abs(entry)) ~= entry | entry == 0
    errordlg('Sample Interval must be a positive integer','Invalid Input','modal')
    set(h,'String',num2str(1e3));
end

% --------------------------------------------------------------------
function varargout = modemenu_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = loadmenu_Callback(h, eventdata, handles, varargin)
% open loadgui and get data
[handles.data,handles.tree] = loadgui(handles.data,handles.tree);
guidata(gcbo,handles);
% update text fields in gui
if ~isempty(handles.data.file)
    set(handles.datafiletxt,'String',strtrunc(handles.data.file,19));
    set(handles.datadirtxt,'String',strtrunc(handles.data.path,19));
    set([handles.numlangtxt,handles.numcogtxt],{'String'},num2cell(size(handles.data.array)'));
end
if  ~isempty(handles.tree.file)
    set(handles.treedirtxt,'String',strtrunc(handles.tree.path,19));
    set(handles.treefiletxt,'String',strtrunc(handles.tree.file,19));
end


% --------------------------------------------------------------------
function varargout = analmenu_Callback(h, eventdata, handles, varargin)

analgui(handles.data,handles.output);
guidata(gcbo,handles);

% --------------------------------------------------------------------
function varargout = specmurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
set([handles.randommurb handles.fixmurb],'Value',0);
set([handles.muvalfixet,handles.muvalet],{'Enable'},{'off';'on'});


% --------------------------------------------------------------------
function varargout = drawtreescb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = plotstatscb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = viewinitbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

if ~isempty(handles.tree.tree)
    draw(handles.tree.tree,handles.output.treefig,LEAF,'Initial Tree');
else
    disp('No valid tree loaded to display');
end

% --------------------------------------------------------------------
function varargout = quitmenu_Callback(h, eventdata, handles, varargin)

delete(gcbf);


% --------------------------------------------------------------------
function varargout = treegui_KeyPressFcn(h, eventdata, handles, varargin)
currchar = get(h,'CurrentChar');
if ~isempty(currchar)
    double(currchar)
end



% --------------------------------------------------------------------
function varargout = savemenu_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = mcoutbutt_Callback(h, eventdata, handles, varargin)

[filename pathname] = uiputfile({'*.nex','Nexus file'},'Select or create an output file');
if ~isequal(filename, 0) & ~isequal(pathname,0)
    if strcmp(strtrunc(filename,4),'.nex') | strcmp(strtrunc(filename,4),'.txt')
        % selected a .nex or .txt. file - remove extension and save
        set([handles.outfiletxt,handles.outdirtxt],{'String'},{filename(1:end-4);pathname});
        handles.output.file = filename(1:end-4);
        handles.output.path = pathname;
    elseif isempty(strfind(filename,'.'))
        % file with no extension - assumes new and overwrite
        set([handles.outfiletxt,handles.outdirtxt],{'String'},{filename;pathname});
        handles.output.file = filename;
        handles.output.path = pathname;
    else
        % not of correct type
        disp('Output file must be of type .txt or .nex')
        disp('No new output file selected')
    end
else
    disp('No new output file selected')
end
guidata(gcbf,handles);

% --------------------------------------------------------------------
function varargout = curr2initbutt_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = viewlastbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

if ~isempty(handles.output.trees)
    pos = length(handles.output.trees);
    s = rnextree(handles.output.trees{pos});
    draw(s,handles.output.treefig,LEAF,['Sample number ' num2str(pos)]);
else
    disp('No output tree to display')
end

% --------------------------------------------------------------------
function varargout = seedrandcb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = quietcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==1
set([handles.drawtreescb handles.plotstatscb],'Enable','off')
else
set([handles.drawtreescb handles.plotstatscb],'Enable','on')
end    



% --------------------------------------------------------------------
function varargout = thetasynet_Callback(h, eventdata, handles, varargin)
entry = str2double(get(h,'string'));
if isnan(entry) | entry <= 0 | entry >= 1
    errordlg('Theta must be a small positive number','Invalid Input','modal')
    set(h,'String',num2str(1/1000));
end


% --------------------------------------------------------------------
function varargout = rootet_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry) | entry <= 0 
    errordlg('Maximum root must be a relatively large positive number','Invalid Input','modal')
    set(h,'String',num2str(160000));
end



% --------------------------------------------------------------------
function varargout = lostonescb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==1
set([handles.lostoneset],'Enable','on')
else
set([handles.lostoneset],'Enable','off')
end    

% --------------------------------------------------------------------
function varargout = lostoneset_Callback(h, eventdata, handles, varargin)
entry = str2double(get(h,'string'));
if isnan(entry) | entry > 1 | entry < 0 
    errordlg('Rare cognates value must be zero or one','Invalid Input','modal')
    set(h,'String',num2str(1));
end

% --------------------------------------------------------------------
function varargout = synthdatacb_Callback(h, eventdata, handles, varargin)

synhand = [handles.langsizeet handles.vocabet handles.synthrandtreerb handles.synthcurrtreerb handles.thetasynet];
    
if get(h,'Value')==0
    set(synhand,'Enable','off');
    if get(handles.issynthcb,'Value')==0
        tt=handles.truetreerb;
        if get(tt,'Value')==1
            set([handles.randtreerb],{'Value'},{1});
            set(handles.initthetaet,'Enable','on');
        end
        set(tt,{'Value','Enable'},{0,'off'});    
        
    end
else
    set([synhand handles.truetreerb],'Enable','on');
    if get(handles.synthrandtreerb,'Value')==0
        set(handles.thetasynet,'Enable','off');
    end
end   

% --------------------------------------------------------------------
function varargout = langsizeet_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry) | floor(entry)~=entry | entry<=0 
    errordlg('Number of languages must be a positive integer','Invalid Input','modal')
    set(h,'String',num2str(40));
end

% --------------------------------------------------------------------
function varargout = vocabet_Callback(h, eventdata, handles, varargin)
entry = str2double(get(h,'string'));
if isnan(entry) | floor(entry)~=entry | entry<=0 
    errordlg('Mean vocabulary size must be a positive integer','Invalid Input','modal')
    set(h,'String',num2str(195));
end

% --------------------------------------------------------------------
function varargout = issynthcb_Callback(h, eventdata, handles, varargin)

tt=handles.truetreerb;
if get(h,'Value')==0
    if get(handles.synthdatacb,'Value')==0
        if get(tt,'Value')==1
            set([handles.randtreerb],{'Value'},{1});
            set(handles.initthetaet,'Enable','on');
        end
        set(tt,{'Value','Enable'},{0,'off'});
    end
else
    set(tt,'Enable','on');
end

% --------------------------------------------------------------------
function varargout = synthrandtreerb_Callback(h, eventdata, handles, varargin)
set(handles.synthcurrtreerb,'Value',0);
set(handles.thetasynet,'Enable','on');

% --------------------------------------------------------------------
function varargout = synthcurrtreerb_Callback(h, eventdata, handles, varargin)
set(handles.thetasynet,'Enable','off');
set(handles.synthrandtreerb,'Value',0);

% --------------------------------------------------------------------
function varargout = truetreerb_Callback(h, eventdata, handles, varargin)

set([handles.spectreerb handles.randtreerb],'Value',0);
set(handles.initthetaet,'Enable','off');

% --------------------------------------------------------------------
function varargout = priorpu_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = initthetaet_Callback(h, eventdata, handles, varargin)
entry = str2double(get(h,'string'));
if isnan(entry) | entry <= 0 | entry >= 1
    errordlg('Theta must be a small positive number','Invalid Input','modal')
    set(h,'String',num2str(1/1000));
end

% --------------------------------------------------------------------
function varargout = syndataoutbutt_Callback(h, eventdata, handles, varargin)

[filename pathname] = uiputfile({'*.nex','Nexus file'},'Select or create an output file');
if ~isequal(filename, 0) & ~isequal(pathname,0)
    if strcmp(strtrunc(filename,4),'.nex')
        % selected a .nex file 
        set([handles.syndataoutfiletxt,handles.syndataoutdirtxt],{'String'},{filename;pathname});
    elseif isempty(strfind(filename,'.'))
        % file with no extension - add extension
        set([handles.syndataoutfiletxt,handles.syndataoutdirtxt],{'String'},{[filename '.nex'];pathname});
    else
        % not of correct type
        disp('Output file must be of type .nex')
        disp('No new output file selected')
    end
else
    disp('No new output file selected')
end
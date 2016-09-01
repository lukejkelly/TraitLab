function varargout = loadgui(varargin)

if isstruct(varargin{1}) % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);    
    handles.data = varargin{1};
    handles.tree = varargin{2};

    guidata(fig, handles);

    % set text fields if already have data files 
    if ~isempty(handles.data.file)
        set(handles.datafiletxt,'String',handles.data.file);
  %      set([handles.langtxt,handles.cogtxt],{'String'},num2cell(size(handles.iodata.array)'));
    end
    if  ~isempty(handles.tree.file)
        set(handles.treefiletxt,'String',handles.tree.file);
    end
         
    % Wait for callbacks to run and window to be dismissed:
    uiwait(fig);
    
    if ~ishandle(fig)
        % Figure dismissed by window control or cancel button.
        % By default, return what we were given
        varargout = varargin;
    else
        % figure dismissed by done button
        % update guidata
        handles = guidata(fig);
        varargout{1} = handles.data;
        varargout{2} = handles.tree;
        delete(fig)
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
function varargout = treecb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
% set([handles.datacb handles.bothcb],'Value',0);

% --------------------------------------------------------------------
function varargout = datacb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
% set([handles.treecb handles.bothcb],'Value',0);

% --------------------------------------------------------------------
% function varargout = bothcb_Callback(h, eventdata, handles, varargin)
% 
% % turn off other radio buttons in group
% set([handles.treecb handles.datacb],'Value',0);

% --------------------------------------------------------------------
function varargout = loadbutt_Callback(h, eventdata, handles, varargin)

readdata = get([handles.datacb],'Value');
readtree = get([handles.treecb],'Value');
usecurr = get([handles.usecurrcb],'Value');

if ~readdata & ~readtree
    errordlg('You must select one of data or tree to load')
    return
end

usecurr = get(handles.usecurrcb,'Value');
if ~usecurr
    % get file names from user
    [filename,pathname] = uigetfile({'*.nex','Nexus'},'Choose a nexus file to load data from');
    if isequal(filename,0)|isequal(pathname,0)
        %no file selected - done 
        disp('No file selected')
        return
    else
        %user selected file - save info in relevant place
        disp(['File ', pathname, filename, ' selected'])
        if readdata
            handles.data.file = filename;
            handles.data.path = pathname;
        end
        if readtree
            handles.tree.file = filename;
            handles.tree.path = pathname;
        end
    end
else
    %check that we have a data or tree file in memory
    if readdata & (isempty(handles.data.file)|isempty(handles.data.path))
        errordlg('There is no current data file loaded.','Error using current files')
        return
    end
        if readtree & (isempty(handles.tree.file)|isempty(handles.tree.path))
        errordlg('There is no current tree file loaded.','Error using current files')
        return
    end
end
    % see if file names are the same
    samefile = strcmp(handles.tree.file,handles.data.file)&strcmp(handles.tree.path,handles.data.path);
if readtree & readdata & samefile
    [s,content] = nexus2stype([handles.tree.path handles.tree.file],readdata,readtree);
    handles.tree.tree = s;
    handles.data.array=content.array;
    handles.data.language=content.language;
    handles.data.cognate=content.cognate;
elseif readtree
    [s,content] = nexus2stype([handles.tree.path handles.tree.file],readdata,readtree);
    handles.tree.tree = s;
elseif readdata
    [s,content] = nexus2stype([handles.data.path handles.data.file],readdata,readtree);
    handles.data.array=content.array;
    handles.data.language=content.language;
    handles.data.cognate=content.cognate;
end
set([handles.treefiletxt,handles.treedirtxt,handles.datafiletxt,handles.datadirtxt],{'String'},{handles.tree.file;handles.tree.path;handles.data.file;handles.data.path});
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = donebutt_Callback(h, eventdata, handles, varargin)

uiresume(gcbf)

% --------------------------------------------------------------------
function varargout = usecurrcb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = cancelbutt_Callback(h, eventdata, handles, varargin)

delete(gcbf)
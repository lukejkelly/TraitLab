function varargout = TraitLab(varargin)
% TraitLab Application M-file for TraitLab.fig

GlobalSwitches;
GlobalValues;
 
if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
        
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    % initialise other data stored in handles
    handles.data  = struct('array',{[]},...
        'language',{''},...
        'cognate',{''},...
        'file',{''},'path',{''},...
        'true',{pop('true')},'truepresent',{0},...
        'clade',{{}},...
        'cladefile',{''},'cladepath',{''});
    handles.tree = struct('tree',{TreeNode([],[],[],[],[],[])},...
        'file',{''},'path',{''},'output',{pop('output')});
    handles.output =pop('output');
    handles.output.file='tloutput';
    if ispc
        handles.output.path=[cd '\'];
    else
        handles.output.path=[cd '/'];
    end
    handles.oldstart.path ='';
    handles.oldstart.file = '';
    guidata(fig, handles);
    
    set([handles.outdirtxt handles.outfiletxt],{'String'},{strtrunc(handles.output.path,19);handles.output.file});
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    
%     try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
%     catch
%         disp(lasterr);
%     end
    
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
able([],[handles.randommurb handles.specmurb],handles.muvalfixet,handles.muvalet);

% --------------------------------------------------------------------
function varargout = randommurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.fixmurb handles.specmurb],[],[handles.muvalfixet,handles.muvalet]);


% --------------------------------------------------------------------
function varargout = muvalet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,0.18);

% --------------------------------------------------------------------
function varargout = muvalfixet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,0.18);

% --------------------------------------------------------------------
function varargout = randtreerb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.spectreerb,handles.truetreerb],[handles.initthetaet handles.maskcb] ,[handles.viewtruebutt]);
if get(handles.maskcb,'Value')
    set(handles.masket,'Enable','on');
end

% --------------------------------------------------------------------
function varargout = spectreerb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.randtreerb,handles.truetreerb],[handles.maskcb],[handles.initthetaet handles.viewtruebutt]);
if get(handles.maskcb,'Value')
    set(handles.masket,'Enable','on');
end
% --------------------------------------------------------------------
function varargout = runet_Callback(h, eventdata, handles, varargin)

checknum(h,1,1e10,1e5,1);


% --------------------------------------------------------------------
function varargout = sampleet_Callback(h, eventdata, handles, varargin)

checknum(h,1,1e6,1e3,1);

% --------------------------------------------------------------------
function varargout = modemenu_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = analmenu_Callback(h, eventdata, handles, varargin)

analgui(handles.data,handles.output);
guidata(gcbo,handles);

% --------------------------------------------------------------------
function varargout = specmurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.randommurb handles.fixmurb],handles.muvalet,handles.muvalfixet);

% --------------------------------------------------------------------
function varargout = drawtreescb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = plotstatscb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = viewinitbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

if handles.tree.output.Nsamp > 0
    pos = str2double(get(handles.numtreeet,'String'));
    s = rnextree(handles.tree.output.trees{pos});
    draw(s,handles.tree.output.treefig,LEAF,sprintf('Tree number %1.0f in file %s',pos,handles.output.file));
else
    disp('No valid tree loaded to display');
end

% --------------------------------------------------------------------
function varargout = quitmenu_Callback(h, eventdata, handles, varargin)

delete(gcbf);

% --------------------------------------------------------------------
function varargout = synthmenu_Callback(h, eventdata, handles, varargin)

synthgui;


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
if get(h,'Value')==0
set([handles.seedet],'Enable','off')
else
set([handles.seedet],'Enable','on')
end    

% --------------------------------------------------------------------
function varargout = quietcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==1
set([handles.drawtreescb handles.plotstatscb],'Enable','off')
else
set([handles.drawtreescb handles.plotstatscb],'Enable','on')
end    

% --------------------------------------------------------------------
function varargout = rootet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1e10,16000);

% --------------------------------------------------------------------
function varargout = lostonescb_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = truetreerb_Callback(h, eventdata, handles, varargin)

%able([],[handles.spectreerb handles.randtreerb],[handles.viewtruebutt],[handles.initthetaet handles.maskcb handles.masket]);
able([],[handles.spectreerb handles.randtreerb],[handles.viewtruebutt],[handles.initthetaet]);

% --------------------------------------------------------------------
function varargout = flatpriorrb_Callback(h, eventdata, handles, varargin)

% turn off yulerb and enable max root age et
able([],handles.yulepriorrb,handles.rootet,[])

% -------------------------------------------------------------------------
function yulepriorrb_Callback(h, eventdata, handles)
% turn off rootrb and disable max root age et
able([],handles.flatpriorrb,[],handles.rootet)

% --------------------------------------------------------------------
function varargout = initthetaet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,1e-3);

% --------------------------------------------------------------------
function varargout = seedet_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry)
    errordlg('Random generator seed must be a real number','Invalid Input','modal')
    set(h,'String',2);
end

% --------------------------------------------------------------------
function varargout = datafilebutt_Callback(h, eventdata, handles, varargin)

% get file names from user
[filename,pathname] = uigetfile({'*.nex','Nexus'},'Choose a nexus file to load data from');
if isequal(filename,0)|isequal(pathname,0)
    %no file selected - done 
    disp('No file selected')
else
    %user selected file - save info in relevant place
    disp(['File ', pathname, filename, ' selected'])
    handles.data.file = filename;
    handles.data.path = pathname;
    set(handles.datafiletxt,'String',strtrunc(handles.data.file,19));
    set(handles.datadirtxt,'String',strtrunc(handles.data.path,19));
    [s,content,handles.data.true,handles.data.clade] = nexus2stype([pathname filename]);
    handles.data.array = content.array;
    handles.data.language = content.language;
    handles.data.cognate = content.cognate;
    handles.data.truepresent = (~isempty(s));
    if handles.data.truepresent
        handles.data.true.state.tree = s;
    else
        handles.data.true.state.tree = [];
    end
        
    [NS,L] = size(handles.data.array);
    set([handles.numlangtxt handles.numcogtxt],{'String'},{NS;L})
    if ~isempty(handles.data.clade)
        able(handles.cladescb,[],handles.cladescb,[]);
        set(handles.cladestxt,'String',sprintf('%1.0f clades found in file',length(handles.data.clade)));
        for hdci=1:length(handles.data.clade)
            disp([sprintf('Clade %3d ',hdci),sprintf(' %s',handles.data.clade{hdci}.language{:})]);
        end
    else
        able([],handles.cladescb,[],handles.cladescb);
        set(handles.cladestxt,'String','');
    end
    if handles.data.truepresent
        set(handles.issynthtxt,'String','File contains synthetic data')
        able([],[],[handles.truetreerb handles.viewtruebutt],[])
    else
        set(handles.issynthtxt,'String','')
        able([],[],[],[handles.truetreerb handles.viewtruebutt])
        if get(handles.truetreerb,'Value')
            able([handles.randtreerb],handles.truetreerb,[handles.maskcb],[]);
            if get(handles.maskcb,'Value')
                set(handles.masket,'Enable','on')
            end
        end
    end           
   guidata(h,handles);
   L=size(handles.data.language,1);
   for k=1:L,
       disp(sprintf('%g %s', k, handles.data.language{k}));
   end
   save outDC;
end



% --------------------------------------------------------------------
function varargout = treefilebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;


[filename pathname] = uigetfile({'*.nex','Nexus'},'Select an output file');

if isequal(filename,0)|isequal(pathname,0)
    %no file selected - done 
    disp('No file selected')
    return
else    
    %user selected file - check that it is a .txt or a .nex file
    if strcmp(strtrunc(filename,4),'.nex') | strcmp(strtrunc(filename,4),'.txt')
        % good file - load output
        disp(sprintf('\nFile %s%s selected',pathname,filename))
        disp('Extracting trees and loss rates')
        [handles.tree.output,ok]=readoutput([pathname, filename(1:end-4)]);
        handles.oldstart.path = pathname;
        handles.oldstart.file = filename;        
        set([handles.treedirtxt handles.treefiletxt],{'String'},{strtrunc(pathname,28);strtrunc(filename,22)})
        guidata(h,handles);
        if ok & handles.tree.output.Nsamp >= 1
            % there is output loaded
            hdl = [handles.numtreeet handles.viewinitbutt];
            set(hdl,'Enable','on');
            set(handles.numtreetxt, 'String',sprintf('The file contains %1.0f trees',handles.tree.output.Nsamp))
            set(handles.numtreeet,'String','1');
            disp(sprintf('%1.0f trees found',handles.tree.output.Nsamp))                
        else
            % problem with loading or no trees in file
            hdl = [handles.numtreeet handles.viewinitbutt];
            set(hdl,'Enable','off');
            disp('No trees loaded')
            set(handles.numtreetxt, 'String','The file contains 0 trees')
         end
    else
        % bad file type
        disp('No file opened - file type must be .nex')
    end
end

% --------------------------------------------------------------------
function varargout = numtreeet_Callback(h, eventdata, handles, varargin)

checknum(h,1,handles.tree.output.Nsamp,1,1);



% --------------------------------------------------------------------
function varargout = maskcb_Callback(h, eventdata, handles, varargin)
if get(h,'Value')
    set(handles.masket,'Enable','on');
else
    set(handles.masket,'Enable','off');
end
% --------------------------------------------------------------------
function varargout = masket_Callback(h, eventdata, handles, varargin)

maskstr = get(h,'String');
if ~isempty(maskstr)
    mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
    if isempty(mask)
        errordlg([maskstr ' is not a valid vector of taxon numbers'],'Invalid Input','modal')
    else
        NS = size(handles.data.array,1);
        if any(floor(mask)~=mask) | mask(1) < 1 | mask(end) > NS
            errordlg(sprintf('Taxa to omit must be a vector of integers between 1 and %1.0f',NS),'Invalid Input','modal')
        end
    end
end

% --------------------------------------------------------------------
function varargout = viewtruebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
if ~isempty(handles.data.true.state.tree)
    draw(handles.data.true.state.tree,handles.output.truefig,LEAF,['True state from ' handles.data.file]);
else
    disp('No true state loaded to display')
end


% --------------------------------------------------------------------
function varargout = cladescb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')
    set(handles.clademasket,'Enable','on');
else
    set(handles.clademasket,'Enable','off');
end

% --------------------------------------------------------------------
function varargout = clademasket_Callback(h, eventdata, handles, varargin)
maskstr = get(h,'String');
if ~isempty(maskstr)
    mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
    if isempty(mask)
        errordlg([maskstr ' is not a valid vector of clade numbers'],'Invalid Input','modal')
    else
        numclade = length(handles.data.clade);
        if any(floor(mask)~=mask) | mask(1) < 1 | mask(end) > numclade
            errordlg(sprintf('Clades to omit must be a vector of integers between 1 and %1.0f',numclade),'Invalid Input','modal')
        end
    end
end

% --------------------------------------------------------------------
function varargout = cogmaskcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')
    set(handles.cogmasket,'Enable','on');
else
    set(handles.cogmasket,'Enable','off');
end

% --------------------------------------------------------------------
function varargout = cogmasket_Callback(h, eventdata, handles, varargin)

maskstr = get(h,'String');
if ~isempty(maskstr)
    mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
    if isempty(mask)
        errordlg([maskstr ' is not a valid vector of trait numbers'],'Invalid Input','modal')
    else
        numcog = size(handles.data.array,2);
        if any(floor(mask)~=mask) | mask(1) < 1 | mask(end) > numcog
            errordlg(sprintf('Traits to omit must be a vector of integers between 1 and %1.0f',numcog),'Invalid Input','modal')
        end
    end
end

% ---------------------------------------------------------------------
function varytopcb_Callback(h, eventdata, handles)


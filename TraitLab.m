function varargout = TraitLab(varargin)
% TraitLab Application M-file for TraitLab.fig

GlobalSwitches;
GlobalValues;
addpath core guifiles borrowing

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

    % Decrease font sizes if we're on a windows computer.
    decreaseFontSizesIfReq(handles)

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
function  pausebutt_Callback(h, eventdata, handles, varargin)

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
function  stopbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
GlobalValues;

set(handles.startbutt,{'Enable','UserData'},{'on',~STOPRUN});
set([handles.pausebutt,h],'Enable','off');
set(handles.pausebutt,'String','Pause');
set(handles.statustxt,'String','Idle');

% --------------------------------------------------------------------
function  fixmurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.randommurb handles.specmurb],handles.muvalfixet,handles.muvalet);

% --------------------------------------------------------------------
function  randommurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.fixmurb handles.specmurb],[],[handles.muvalfixet,handles.muvalet]);


% --------------------------------------------------------------------
function  muvalet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,0.18);

% --------------------------------------------------------------------
function  muvalfixet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,0.18);

% --------------------------------------------------------------------
function  randtreerb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.spectreerb,handles.truetreerb],[handles.initthetaet handles.maskcb] ,[handles.viewtruebutt]);
if get(handles.maskcb,'Value')
    set(handles.masket,'Enable','on');
end

% --------------------------------------------------------------------
function  spectreerb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.randtreerb,handles.truetreerb],[handles.maskcb],[handles.initthetaet handles.viewtruebutt]);
if get(handles.maskcb,'Value')
    set(handles.masket,'Enable','on');
end
% --------------------------------------------------------------------
function  runet_Callback(h, eventdata, handles, varargin)

checknum(h,1,1e10,1e5,1);


% --------------------------------------------------------------------
function sampleet_Callback(h, eventdata, handles, varargin)

checknum(h,1,1e6,1e3,1);

% --------------------------------------------------------------------
function  modemenu_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function  analmenu_Callback(h, eventdata, handles, varargin)

analgui(handles.data,handles.output);
guidata(gcbo,handles);

% --------------------------------------------------------------------
function specmurb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.randommurb handles.fixmurb],handles.muvalet,handles.muvalfixet);

% --------------------------------------------------------------------
function  drawtreescb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function  plotstatscb_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function  viewinitbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

if handles.tree.output.Nsamp > 0
    pos = str2double(get(handles.numtreeet,'String'));
    s = rnextree(handles.tree.output.trees{pos});
    %keyboard;
    if length(handles.tree.output.cattrees)==length(handles.tree.output.trees)
        %they should see that it is a catree even if they dont plan to use
        %the cats - GKN 18/3/11
        sc=rnextree(handles.tree.output.cattrees{pos});
        %[s, catlist]=CatTreeToList(sc,s);
		s = CatTreeToList(sc,s);
    end
    draw(s,handles.tree.output.treefig,LEAF,sprintf('Tree number %1.0f in file %s',pos,handles.output.file));
else
    disp('No valid tree loaded to display');
end

% --------------------------------------------------------------------
function  quitmenu_Callback(h, eventdata, handles, varargin)

delete(gcbf);

% --------------------------------------------------------------------
function  synthmenu_Callback(h, eventdata, handles, varargin)

synthgui;


% --------------------------------------------------------------------
function  mcoutbutt_Callback(h, eventdata, handles, varargin)

[filename pathname] = uiputfile({'*.nex','Nexus file'},'Select or create an output file');
if ~isequal(filename, 0) && ~isequal(pathname,0)
    if strcmp(strtrunc(filename,4),'.nex') || strcmp(strtrunc(filename,4),'.txt')
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
function viewlastbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

if ~isempty(handles.output.trees)
    pos = length(handles.output.trees);
    s = rnextree(handles.output.trees{pos});
    if MCMCCAT
        scat= rnextree(handles.output.cattrees{pos});
        s=CatTreeToList(scat,s);
    end
    draw(s,handles.output.treefig,LEAF,['Sample number ' num2str(pos)]);
else
    disp('No output tree to display')
end

% --------------------------------------------------------------------
function  seedrandcb_Callback(h, eventdata, handles, varargin)
if get(h,'Value')==0
set([handles.seedet],'Enable','off')
else
set([handles.seedet],'Enable','on')
end

% --------------------------------------------------------------------
function quietcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==1
set([handles.drawtreescb handles.plotstatscb],'Enable','off')
else
set([handles.drawtreescb handles.plotstatscb],'Enable','on')
end

% --------------------------------------------------------------------
function  rootet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1e10,16000);

% --------------------------------------------------------------------
function  lostonescb_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function truetreerb_Callback(h, eventdata, handles, varargin)

%able([],[handles.spectreerb handles.randtreerb],[handles.viewtruebutt],[handles.initthetaet handles.maskcb handles.masket]);
able([],[handles.spectreerb handles.randtreerb],[handles.viewtruebutt],[handles.initthetaet]);

% --------------------------------------------------------------------
function   flatpriorrb_Callback(h, eventdata, handles, varargin)

% turn off yulerb and enable max root age et
able([],handles.yulepriorrb,handles.rootet,[])

% -------------------------------------------------------------------------
function yulepriorrb_Callback(h, eventdata, handles)
% turn off rootrb and disable max root age et
able([],handles.flatpriorrb,[],handles.rootet)

% --------------------------------------------------------------------
function   initthetaet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,1e-3);

% --------------------------------------------------------------------
function  seedet_Callback(h, eventdata, handles, varargin)

entry = str2double(get(h,'string'));
if isnan(entry)
    errordlg('Random generator seed must be a real number','Invalid Input','modal')
    set(h,'String',2);
end

% --------------------------------------------------------------------
function  datafilebutt_Callback(h, eventdata, handles, varargin)

% get file names from user
[filename,pathname] = uigetfile({'*.nex','Nexus'},'Choose a nexus file to load data from');
if isequal(filename,0)||isequal(pathname,0)
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
       fprintf('%g %s\n', k, handles.data.language{k});
   end
   %% DW 19/7/2007
   % report those languages missing at least 5% of data
   % [NS,L] = size(handles.data.array);
   %missingmany = find(sum(handles.data.array == 2,2) >= L/20);
   missingmany = find(sum(handles.data.array == 2,2) >= 1);
   %length(missingmany)
   if ~isempty(missingmany)
       fprintf('\nThe following taxa have at least 5%% missing data: \n')
       for k = missingmany'
           %disp(sprintf('%g %s %d%%', k, handles.data.language{k},round(100*sum(handles.data.array(k,:) == 2)/L)))
           fprintf('%g %s %d\n', k, handles.data.language{k},round(sum(handles.data.array(k,:) == 2)))
       end
   end
   %% DW 19/7/2007 end


   save outDC;
end



% --------------------------------------------------------------------
function  treefilebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;


[filename pathname] = uigetfile({'*.nex','Nexus'},'Select an output file');

if isequal(filename,0)||isequal(pathname,0)
    %no file selected - done
    disp('No file selected')
    return
else
    %user selected file - check that it is a .txt or a .nex file
    if strcmp(strtrunc(filename,4),'.nex') || strcmp(strtrunc(filename,4),'.txt')
        % good file - load output
        fprintf('\nFile %s%s selected\n',pathname,filename)
        disp('Extracting trees and loss rates')
        [handles.tree.output,ok]=readoutput([pathname, filename(1:end-4)]);
        handles.oldstart.path = pathname;
        handles.oldstart.file = filename;
        set([handles.treedirtxt handles.treefiletxt],{'String'},{strtrunc(pathname,28);strtrunc(filename,22)})
        guidata(h,handles);
        if ok && handles.tree.output.Nsamp >= 1
            % there is output loaded
            hdl = [handles.numtreeet handles.viewinitbutt];
            set(hdl,'Enable','on');
            set(handles.numtreetxt, 'String',sprintf('The file contains %1.0f trees',handles.tree.output.Nsamp))
            set(handles.numtreeet,'String','1');
            fprintf('%1.0f trees found\n',handles.tree.output.Nsamp)
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
function numtreeet_Callback(h, eventdata, handles, varargin)

checknum(h,1,handles.tree.output.Nsamp,1,1);



% --------------------------------------------------------------------
function  maskcb_Callback(h, eventdata, handles, varargin)
if get(h,'Value')
    set(handles.masket,'Enable','on');
else
    set(handles.masket,'Enable','off');
end
% --------------------------------------------------------------------
function  masket_Callback(h, eventdata, handles, varargin)

maskstr = get(h,'String');
if ~isempty(maskstr)
    mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
    if isempty(mask)
        errordlg([maskstr ' is not a valid vector of taxon numbers'],'Invalid Input','modal')
    else
        NS = size(handles.data.array,1);
        if any(floor(mask)~=mask) || mask(1) < 1 || mask(end) > NS
            errordlg(sprintf('Taxa to omit must be a vector of integers between 1 and %1.0f',NS),'Invalid Input','modal')
        end
    end
end

% --------------------------------------------------------------------
function  viewtruebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;
if ~isempty(handles.data.true.state.tree)
    draw(handles.data.true.state.tree,handles.output.truefig,LEAF,['True state from ' handles.data.file]);
else
    disp('No true state loaded to display')
end


% --------------------------------------------------------------------
function  cladescb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')
    set(handles.clademasket,'Enable','on');
else
    set(handles.clademasket,'Enable','off');
end

% --------------------------------------------------------------------
function  clademasket_Callback(h, eventdata, handles, varargin)
maskstr = get(h,'String');
if ~isempty(maskstr)
    mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
    if isempty(mask)
        errordlg([maskstr ' is not a valid vector of clade numbers'],'Invalid Input','modal')
    else
        numclade = length(handles.data.clade);
        if any(floor(mask)~=mask) || mask(1) < 1 || mask(end) > numclade
            errordlg(sprintf('Clades to omit must be a vector of integers between 1 and %1.0f',numclade),'Invalid Input','modal')
        end
    end
end

% --------------------------------------------------------------------
function  cogmaskcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')
    set(handles.cogmasket,'Enable','on');
else
    set(handles.cogmasket,'Enable','off');
end

% --------------------------------------------------------------------
function  cogmasket_Callback(h, eventdata, handles, varargin)

maskstr = get(h,'String');
if ~isempty(maskstr)
    mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
    if isempty(mask)
        errordlg([maskstr ' is not a valid vector of trait numbers'],'Invalid Input','modal')
    else
        numcog = size(handles.data.array,2);
        if any(floor(mask)~=mask) || mask(1) < 1 || mask(end) > numcog
            errordlg(sprintf('Traits to omit must be a vector of integers between 1 and %1.0f',numcog),'Invalid Input','modal')
        end
    end
end

% ---------------------------------------------------------------------
function varytopcb_Callback(h, eventdata, handles)



% --- Executes on button press in fixkapparb.
function fixkapparb_Callback(hObject, eventdata, handles)
% turn off other radio buttons in group
able([],[handles.randomkapparb handles.speckapparb],handles.kappavalfixet,handles.kappavalet);


% --- Executes on button press in randomkapparb.
function randomkapparb_Callback(hObject, eventdata, handles)
% turn off other radio buttons in group
able([],[handles.fixkapparb handles.speckapparb],[],[handles.kappavalfixet handles.kappavalet]);



function kappavalet_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function kappavalet_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in speckapparb.
function speckapparb_Callback(hObject, eventdata, handles)
% turn off other radio buttons in group
able([],[handles.randomkapparb handles.fixkapparb],handles.kappavalet,handles.kappavalfixet);



function kappavalfixet_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of kappavalfixet as text
%        str2double(get(hObject,'String')) returns contents of kappavalfixet as a double


% --- Executes during object creation, after setting all properties.
function kappavalfixet_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fixrhorb.
function fixrhorb_Callback(hObject, eventdata, handles)
% turn off other radio buttons in group
able([],[handles.randomrhorb handles.specrhorb],handles.rhovalfixet,handles.rhovalet);


% --- Executes on button press in randomrhorb.
function randomrhorb_Callback(hObject, eventdata, handles)
% turn off other radio buttons in group
able([],[handles.fixrhorb handles.specrhorb],[],[handles.rhovalfixet handles.rhovalet]);



function rhovalet_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of rhovalet as text
%        str2double(get(hObject,'String')) returns contents of rhovalet as a double


% --- Executes during object creation, after setting all properties.
function rhovalet_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in specrhorb.
function specrhorb_Callback(hObject, eventdata, handles)
% turn off other radio buttons in group
able([],[handles.randomrhorb handles.fixrhorb],handles.rhovalet,handles.rhovalfixet);



function rhovalfixet_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of rhovalfixet as text
%        str2double(get(hObject,'String')) returns contents of rhovalfixet as a double


% --- Executes during object creation, after setting all properties.
function rhovalfixet_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in includecatscb.
function includecatscb_Callback(h, eventdata, handles)

if get(h,'Value')==0
	%% disable all catastrophe options
	set([handles.fixkapparb,handles.kappavalfixet,handles.speckapparb,handles.kappavalet,handles.randomkapparb],'Enable','off')
	set([handles.fixrhorb,handles.rhovalfixet,handles.specrhorb,handles.rhovalet,handles.randomrhorb],'Enable','off')
else
	%% enable all catastrophe radio buttons
	set([handles.fixkapparb,handles.speckapparb,handles.randomkapparb],'Enable','on')
	set([handles.fixrhorb,handles.specrhorb,handles.randomrhorb],'Enable','on')
	%% enable edit text only with currently selected radio button
	if get(handles.fixkapparb,'Value')==1
		set(handles.kappavalfixet,'Enable','on')
	end
	if get(handles.speckapparb,'Value')==1
		set(handles.kappavalet,'Enable','on')
	end
	if get(handles.fixrhorb,'Value')==1
		set(handles.rhovalfixet,'Enable','on')
	end
	if get(handles.specrhorb,'Value')==1
		set(handles.rhovalet,'Enable','on')
	end
end







% --- Executes on button press in cladeagescb.
function cladeagescb_Callback(h, eventdata, handles)
if get(h,'Value')
    set(handles.cladeagesmasket,'Enable','on');
else
    set(handles.cladeagesmasket,'Enable','off');
end



function cladeagesmasket_Callback(h, eventdata, handles)
maskstr = get(h,'String');
if ~isempty(maskstr)
    mask = sort(unique(str2num(maskstr))); %#ok<ST2NM>
    if isempty(mask)
        errordlg([maskstr ' is not a valid vector of clade numbers'],'Invalid Input','modal')
    else
        numclade = length(handles.data.clade);
        if any(floor(mask)~=mask) || mask(1) < 1 || mask(end) > numclade
            errordlg(sprintf('Clades to omit must be a vector of integers between 1 and %1.0f',numclade),'Invalid Input','modal')
        end
    end
end





% --- Executes during object creation, after setting all properties.
function clademasket_CreateFcn(h, eventdata, handles)
% hObject    handle to clademasket (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(h,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(h,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function cladeagesmasket_CreateFcn(h, eventdata, handles)
% hObject    handle to cladeagesmasket (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(h,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(h,'BackgroundColor','white');
end


% --- Executes on button press in allowForLateralTransferCB.
function allowForLateralTransferCB_Callback(hObject, eventdata, handles)
% hObject    handle to allowForLateralTransferCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allowForLateralTransferCB

if get(hObject, 'Value')
   able([], [], [handles.varyLateralTransferRateCB, ...
       handles.initialiseLateralTransferRateAtSpecifiedValueRB, ...
       handles.initialiseLateralTransferRateAtRandomValueRB, ...
       handles.initialiseLateralTransferRateAtSpecifiedValueEB], []);
else
   able([], [], [], [handles.varyLateralTransferRateCB, ...
       handles.initialiseLateralTransferRateAtSpecifiedValueRB, ...
       handles.initialiseLateralTransferRateAtRandomValueRB, ...
       handles.initialiseLateralTransferRateAtSpecifiedValueEB]);
end


% --- Executes on button press in varyLateralTransferRateCB.
function varyLateralTransferRateCB_Callback(hObject, eventdata, handles)
% hObject    handle to varyLateralTransferRateCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of varyLateralTransferRateCB


% --- Executes on button press in initialiseLateralTransferRateAtSpecifiedValueRB.
function initialiseLateralTransferRateAtSpecifiedValueRB_Callback(hObject, eventdata, handles)
% hObject    handle to initialiseLateralTransferRateAtSpecifiedValueRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of initialiseLateralTransferRateAtSpecifiedValueRB
if get(hObject, 'Value')
    able([], [], handles.initialiseLateralTransferRateAtSpecifiedValueEB, []);
end

% --- Executes on button press in initialiseLateralTransferRateAtRandomValueRB.
function initialiseLateralTransferRateAtRandomValueRB_Callback(hObject, eventdata, handles)
% hObject    handle to initialiseLateralTransferRateAtRandomValueRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of initialiseLateralTransferRateAtRandomValueRB

if get(hObject, 'Value')
    able([], [], [], handles.initialiseLateralTransferRateAtSpecifiedValueEB);
end


function initialiseLateralTransferRateAtSpecifiedValueEB_Callback(hObject, eventdata, handles)
% hObject    handle to initialiseLateralTransferRateAtSpecifiedValueEB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initialiseLateralTransferRateAtSpecifiedValueEB as text
%        str2double(get(hObject,'String')) returns contents of initialiseLateralTransferRateAtSpecifiedValueEB as a double

checknum(hObject, 0, 1, 1e-3);




% --- Executes during object creation, after setting all properties.
function initialiseLateralTransferRateAtSpecifiedValueEB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initialiseLateralTransferRateAtSpecifiedValueEB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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

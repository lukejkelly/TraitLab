function varargout = synthgui(varargin)
% SYNTHGUI Application M-file for synthgui.fig
%    FIG = SYNTHGUI launch synthgui GUI.
%    SYNTHGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 23-May-2011 16:01:38
GlobalSwitches;
GlobalValues;
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    fig.Resize = 'on';  % LJK 12/3/20

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it.
	handles = guihandles(fig);
    % store any loaded trees as a cell array of tree strings
    handles.tree = {};
    handles.cattree={};
    handles.catsmissing=0;
    handles.numtrees = 0;
    handles.path = [cd '/'];
    handles.file = 'synthdata.nex';
	guidata(fig, handles);

    set([handles.datafiletxt handles.datadirtxt],{'String'},{handles.file;handles.path})

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


% --------------------------------------------------------------------
function thetaet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,0.001);

% --------------------------------------------------------------------
function langet_Callback(h, eventdata, handles, varargin)

checknum(h,2,200000,50,1);
% --------------------------------------------------------------------
function vocabet_Callback(h, eventdata, handles, varargin)

checknum(h,0,200000,200,1);

% --------------------------------------------------------------------
function randrb_Callback(h, eventdata, handles, varargin)

able([],[handles.filerb],[handles.thetaet handles.langet],[handles.treefilebutt handles.numtreeet handles.viewtreebutt]);

% --------------------------------------------------------------------
function filerb_Callback(h, eventdata, handles, varargin)

if handles.numtrees > 0
able([],[handles.randrb],[handles.treefilebutt handles.numtreeet handles.viewtreebutt],[handles.thetaet handles.langet]);
else
  able([],[handles.randrb],[handles.treefilebutt],[handles.thetaet handles.langet]);
end

% --------------------------------------------------------------------
function datafilebutt_Callback(h, eventdata, handles, varargin)

[filename pathname] = uiputfile({'*.nex','Nexus file'},'Select or create an output file');
if ~isequal(filename, 0) && ~isequal(pathname,0)
    if strcmp(strtrunc(filename,4),'.nex')
        % selected a .nex file
        set([handles.datafiletxt,handles.datadirtxt],{'String'},{filename;pathname});
        handles.path = pathname;
        handles.file = filename;
    elseif isempty(strfind(filename,'.'))
        % file with no extension - add extension
        set([handles.datafiletxt,handles.datadirtxt],{'String'},{[filename '.nex'];pathname});
                handles.path = pathname;
        handles.file = [filename '.nex'];

    else
        % not of correct type
        disp('Output file must be of type .nex')
        disp('No new output file selected')
    end
else
    disp('No new output file selected')
end

guidata(h,handles);
% --------------------------------------------------------------------
function muet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1,0.18);

% --------------------------------------------------------------------
function treefilebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

[filename,pathname] = uigetfile({'*.nex','Nexus'},'Choose a nexus file containing trees');
if isequal(filename,0)||isequal(pathname,0)
    %no file selected
    disp('No file selected')
else
    %user selected file - save info in relevant place
    fprintf('\nFile %s%s selected.\n',pathname,filename)
    disp('Extracting trees from file')
    handles.tree = readalltrees([pathname, filename]);
    handles.numtrees = length(handles.tree);
    if exist([pathname, filename(1:end-4), 'cat.nex'],'file')
        handles.cattree = readalltrees([pathname, filename(1:end-4), 'cat.nex']);
        if length(handles.cattree)~=handles.numtrees
            warning('Number of trees is different than number on catastrophe trees. This may lead to errors.')
        end
    else
        handles.catsmissing=1;
    end
    guidata(h,handles);
    fprintf('%1.0f trees found in file\n',handles.numtrees);
    set([handles.treefiletxt,handles.treedirtxt],{'String'},{filename;pathname});
    set(handles.numtreetxt,'String',sprintf('The file contains %1.0f trees',handles.numtrees));
    if handles.numtrees >0
        set([handles.numtreeet handles.viewtreebutt],'Enable','on');
        set(handles.numtreeet,'String','1')
    else
                set([handles.numtreeet handles.viewtreebutt],'Enable','off');
    end
end

% --------------------------------------------------------------------
function numtreeet_Callback(h, eventdata, handles, varargin)

checknum(h,1,handles.numtrees,1,1);

% --------------------------------------------------------------------
function borrowet_Callback(h, eventdata, handles, varargin)

checknum(h,0,inf,0.15);

% --------------------------------------------------------------------
function borrowcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==0
    set([handles.borrowet,handles.borrowtxt,handles.localborrowcb,handles.localborrowet,handles.localborrowtxt],'Enable','off')
else
    if get(handles.localborrowcb,'Value')==1
        set([handles.borrowet,handles.borrowtxt,handles.localborrowcb,handles.localborrowet,handles.localborrowtxt],'Enable','on')
    else
        set([handles.borrowet,handles.localborrowcb,handles.borrowtxt],'Enable','on')
    end
end
% --------------------------------------------------------------------
function localborrowet_Callback(h, eventdata, handles, varargin)

checknum(h,0,1e8,1e3);


% --------------------------------------------------------------------
function localborrowcb_Callback(h, eventdata, handles, varargin)
if get(h,'Value')==0
set([handles.localborrowet,handles.localborrowtxt],'Enable','off')
else
set([handles.localborrowet,handles.localborrowtxt],'Enable','on')
end

% --------------------------------------------------------------------
function polymorphet_Callback(h, eventdata, handles, varargin)

checknum(h,1,1e6,180);

% --------------------------------------------------------------------
function polymorphcb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==0
    able([],[],[],[handles.polymorphet,handles.polymorphtxt,handles.rhclasscb,handles.rhclasset,handles.rhclasstxt]);
else
    if get(handles.rhclasscb,'Value')==1
        able([],[],[handles.polymorphet,handles.polymorphtxt,handles.rhclasscb,handles.rhclasset,handles.rhclasstxt],[]);
    else
        able([],[],[handles.polymorphet,handles.polymorphtxt,handles.rhclasscb],[]);
    end
end


% --------------------------------------------------------------------
function viewtreebutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

pos = str2double(get(handles.numtreeet,'string'));
filename = get([handles.treedirtxt, handles.treefiletxt],'String');
filename = [filename{:}];
s = rnextree(handles.tree{pos});
if get(handles.includecatscb,'Value') && ~handles.catsmissing
    MCMCCAT=1;
    scat= rnextree(handles.cattree{pos});
    s=CatTreeToList(scat,s);
else
    MCMCCAT=0;
end
fig = pop('output');
draw(s,fig.treefig,LEAF,sprintf('Tree number %1.0f in %s',pos,filename));
% --------------------------------------------------------------------
function lostcb_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function numcladeset_Callback(h, eventdata, handles, varargin)

checknum(h,1,inf,1,1)



% --------------------------------------------------------------------
function cladeboundet_Callback(h, eventdata, handles, varargin)

checknum(h,0,100,10,0)


% --------------------------------------------------------------------
function cladescb_Callback(h, eventdata, handles, varargin)

if get(h,'Value')==0
able([],[],[],[handles.originrb,handles.cladeboundet,handles.numcladeset,handles.rootrb,handles.bothrb]);
else
able([],[],[handles.originrb,handles.cladeboundet,handles.numcladeset,handles.rootrb,handles.bothrb],[]);
end



% --------------------------------------------------------------------
function originrb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.rootrb handles.bothrb],[],[]);




% --------------------------------------------------------------------
function rootrb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.originrb handles.bothrb],[],[]);



% --------------------------------------------------------------------
function bothrb_Callback(h, eventdata, handles, varargin)
% turn off other radio buttons in group
able([],[handles.originrb handles.rootrb],[],[]);

% --------------------------------------------------------------------

function rhbranchcb_Callback(h, eventdata, handles)

if get(h,'Value')==0
    able([],[],[],[handles.rhbranchet,handles.rhbranchtxt]);
else
    able([],[],[handles.rhbranchet,handles.rhbranchtxt],[]);
end


% --------------------------------------------------------------------

function rhbranchet_Callback(h, eventdata, handles)
% hObject    handle to rhbranchet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhbranchet as text
%        str2double(get(hObject,'String')) returns contents of rhbranchet as a double

checknum(h,0,inf,1,0)

% --------------------------------------------------------------------

% --- Executes on button press in rhclasscb.
function rhclasscb_Callback(h, eventdata, handles)
% hObject    handle to rhclasscb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rhclasscb

if get(h,'Value')==0
    able([],[],[],[handles.rhclasset,handles.rhclasstxt]);
else
    able([],[],[handles.rhclasset,handles.rhclasstxt],[]);
end



% --------------------------------------------------------------------

function rhclasset_Callback(h, eventdata, handles)
% hObject    handle to rhclasset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhclasset as text
%        str2double(get(hObject,'String')) returns contents of rhclasset as
%        a double

checknum(h,0,inf,1,0)

% --------------------------------------------------------------------

function includecatscb_Callback(h, eventdata, handles)

if get(h,'Value')==0
    able([],[],[],[handles.kappaet,handles.rhoet]);
else
    able([],[],[handles.kappaet,handles.rhoet],[]);
end

function rhoet_Callback(h, eventdata, handles)
% hObject    handle to rhclasset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhclasset as text
%        str2double(get(hObject,'String')) returns contents of rhclasset as
%        a double

checknum(h,0,inf,1e-4,0)

function kappaet_Callback(h, eventdata, handles)
% hObject    handle to rhclasset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rhclasset as text
%        str2double(get(hObject,'String')) returns contents of rhclasset as
%        a double

checknum(h,0,1,0.5,0)




function synthbutt_Callback(h, eventdata, handles, varargin)

GlobalSwitches;

% get user values from GUI
NS = 0;       %NUMSEQ number of leaves
TH = 0;       %THETA language branching rate for synthetic data
MK = OFF;     %MASKING switch {ON OFF}
LO = OFF;     %LOSTONES cognates present in just one language are not observed
if get(handles.lostcb,'Value')
    LT = 1;   %delete cognates present in LOST or less of the languages
else
    LT = 0;
end
LR = str2double(get(handles.muet,'String'));  %LOSSRATE proportion cognates lost in 1000 years
BF =  0;   % BORROWFRAC borrowing rate as a fraction of the death rate mu
LB = OFF;  % LOCALBORROW ON for local borrowing (also need maxdist) or OFF for wide borrowing
MD = 0;    % MAXDISTance that languages can borrow within
if get(handles.borrowcb,'Value')
    BO = ON;   %BORROW switch (ON, OFF) simulate borrowing
    BF =  str2double(get(handles.borrowet,'String'));
    if get(handles.localborrowcb,'Value')
        LB = ON;
        MD = str2double(get(handles.localborrowet,'String'));
    end
else
    BO = OFF;
end
if get(handles.polymorphcb,'Value')
    % using polymorphic model - get value for number of meaning classes
    NMC =  str2double(get(handles.polymorphet,'String'));
    PM = ON;
    % see whether there is rate heterogeneity across classes
    if get(handles.rhclasscb,'Value')
        % rate heterogeneity across classes - get variance parameter
        RHC =  str2double(get(handles.rhclasset,'String'));
    else
        RHC =  0;
    end
else
    NMC = 1;   % standard model
    PM = OFF;
    RHC = 0;
end

if get(handles.missingdatacb,'Value')
    MISDAT=1;
    % MCMCMISS=1;
else
    MISDAT=0;
    % MCMCMISS=0;
end

if get(handles.includecatscb,'Value')
    MCMCCAT=1;
%     if handles.catsmissing && get(handles.filerb)
%         disp('Warning: There is no catastrophe file associated with the tree file you selected. Catastrophes will not be included.')
%         MCMCCAT=0;
%     end
else
    MCMCCAT=0;
end

if get(handles.rhbranchcb,'Value')
    % rate heterogeneity on branches - get variance parameter
    RHB =  str2double(get(handles.rhbranchet,'String'));
else
    RHB =  0;
end

if MCMCCAT
    RHO=str2double(get(handles.rhoet,'String'));
    KAPPA=str2double(get(handles.kappaet,'String'));
else
    RHO=0;
    KAPPA=0;
end

VS = str2double(get(handles.vocabet,'String'));  %VOCABSIZE mean number of words in a language
ST  = NEWTRE;    %SYNTHSTYLE {OLDTRE NEWTRE} if DATASOURCE=BUILD
STF = '';        %SYNTHTREFILE name, if SYNTHSTYLE=OLDTRE
STR = '';        %SYNTHTREE
CK=0; % Catastrophes are not already known
numtree = 0;
if get(handles.filerb,'Value')
    ST = OLDTRE;
    STF = get([handles.treedirtxt, handles.treefiletxt],'String');
    STF = [STF{:}];
    if  strmatch(STF,'None selected')
        disp('Need to specify a file with trees in it')
        return
    end
    numtree = str2double(get(handles.numtreeet,'String'));
    STR = rnextree(handles.tree{numtree});
    if get(handles.includecatscb,'Value')
        if handles.catsmissing
            disp('No catastrophe file associated with tree file you selected. New catastrophes will be simulated.')
        else
            STR=CatTreeToList(rnextree(handles.cattree{numtree}),STR);
            CK=1;
        end
    end


else
    NS = str2double(get(handles.langet,'String'));
    TH = str2double(get(handles.thetaet,'String'));
end

%write the control variables into structures used by fullsetup
fsu=pop('fsu');
fsu.RUNLENGTH         = 1     ;
fsu.SUBSAMPLE         = 1     ;
fsu.SEEDRAND          = OFF   ;
fsu.SEED              = 2     ;
fsu.DATASOURCE        = BUILD ;
fsu.DATAFILE          = ''    ;
fsu.DATASYN           = OFF   ;
fsu.SYNTHSTYLE        = ST    ;
fsu.SYNTHTREFILE      = STF   ;
fsu.SYNTHTRE          = STR   ;
fsu.KNOWCATS          = CK    ;
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
fsu.LOSTONES          = LO    ;
fsu.LOST              = LT    ;
fsu.LOSSRATE          = LR    ;
fsu.LOSSRATEBRANCHVAR = RHB   ;
fsu.LOSSRATECLASSVAR  = RHC   ;
fsu.PSURVIVE          = 1     ;
fsu.BORROW            = BO    ;
fsu.BORROWFRAC        = BF    ;
fsu.POLYMORPH         = PM    ;
fsu.LOCALBORROW       = LB    ;
fsu.MAXDIST           = MD    ;
fsu.NMEANINGCLASS     = NMC   ;
fsu.MASKING           = MK    ;
fsu.DATAMASK          = []    ;
fsu.ISCOLMASK         = OFF   ;
fsu.COLUMNMASK        = []    ;
fsu.NUMSEQ            = NS    ;
fsu.VOCABSIZE         = VS    ;
fsu.THETA             = TH    ;
fsu.ISCLADE           = OFF   ;
fsu.CLADE             = []    ;
fsu.STRONGCLADES      = []    ;
fsu.GUITRUE           = []    ;
fsu.GUICONTENT        = []    ;
fsu.SYNTHMISS         =MISDAT ;
fsu.MCMCCAT           =MCMCCAT;
fsu.RHO               =RHO;
fsu.KAPPA             =KAPPA;
fsu.MISDAT            =MISDAT;

data = fullsetup(fsu);


% synthesize clades
clades = {};
if get(handles.cladescb,'Value')
    numclade = str2double(get(handles.numcladeset,'String'));
    % check that the number of clades is less than the number of internal nodes
    if numclade > (data.true.NS - 2)
        fprintf('Number of clades must be less than or equal to number of internal nodes (%1.0f).\n No clades generated\n',data.true.NS-1);
    else
        fprintf('Generating %1.0f clades\n',numclade);
        rbvals=[handles.originrb handles.rootrb handles.bothrb];
        rbvals = get(rbvals,'Value');
        originrootboth =  find([rbvals{:}]);
        accuracy = str2double(get(handles.cladeboundet,'String'))/100;
        clades = synthclades(data.true.state.tree,numclade,originrootboth,accuracy);
        % display clades on screen
        sclde = Clade2Tree(clades);
        fig = pop('output');
        % draw(sclde,fig.cladefig,CONC,'Clade structure schematic'); Not useful: the clades are shown on the complete tree. RJR 23/05/11

        for i = 1:length(clades)
            fprintf('Clade %g:',i)
            disp(clades{i})
        end
    end
end

% save it
fprintf('Saving data to %s%s\n',handles.path, handles.file);
[str,ok] = stype2nexus(data.true.state.tree,['Synthetic data made from Synthesize GUI on ' datestr(clock)],'TRUE',data.true,clades);
if ok
    fid = fopen([handles.path, handles.file],'w');
    if fid > 0
        fprintf(fid,str);
        if fclose(fid) ~= 0
            fprintf('\nProblem with closing %s%s\n',handles.path, handles.file);
        end
        % write details in synthsized data in par file
        parfilename = strrep([handles.path, handles.file],'.nex','.par');
        fprintf('Writing parameter file %s\n',parfilename);
        ok = writesynparfile(parfilename,fsu,numtree);
        if ~ok
            fprintf('Encountered problem writing parameter file %s.\n Data file should be ok.\n',parfilename);
        end
    else
        fprintf('\nProblem with opening %s%s\n  Data not saved\n',handles.path, handles.file);
    end
else
    fprintf('\nSynthetic data not saved\n');
end
fprintf('Data synthesizing complete\n')

fig = pop('output');
draw(data.true.state.tree,fig.truefig,0,'Tree used for synthesising data',[],[],clades);

%[x,y]=DepthDist(data.content,data.true.state);
%figure;plot(x,y,'.');xy=axis; hold on; plot([xy(1),xy(2)],[xy(1),xy(2)]); hold off;

 if get(handles.exploredatacb,'Value')
    exploredata(data.true.state,data.content,1,0,1); %GKN 7/7 dropped useless filename and HEt args
 end


% --- Executes on button press in exploredatacb.
function exploredatacb_Callback(hObject, eventdata, handles)
% hObject    handle to exploredatacb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exploredatacb

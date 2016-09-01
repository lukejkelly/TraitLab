function [makeHTML HTMLname]=gethtmlvar(handles)

makeHTML=get(handles.makeHTML,'Value');

HTMLname=get(handles.HTMLname,'Value');

if makeHTML && ~HTMLname
    % use default name
    HTMLname=datestr(now,'yyyymmddHHMMSS');
end

if makeHTML && strcmpi(HTMLname(end-4:end),'.html')
    % remove .html extension (to avoid it appearing twice)
    HTMLname=HTMLname(1:end-5);
end
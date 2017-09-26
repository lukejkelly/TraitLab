function able(switchon,switchoff,enable,disable)

% able(switchon,switchoff,enable,disable);
%
% All variables passed in are vectors of handles.
% switchon/switchoff have 'Value' set to 1/0,
% enable/disable have 'Enable' set to 'on'/'off'.

if ~isempty(switchon)
    set(switchon,'Value',1);
end

if ~isempty(switchoff)
    set(switchoff,'Value',0);
end

if ~isempty(enable)
    set(enable,'Enable','on');
end

if ~isempty(disable)
    set(disable,'Enable','off');
end
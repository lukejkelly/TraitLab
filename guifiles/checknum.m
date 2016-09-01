function checknum(h,lowerlim,upperlim,default,isint)

% checknum(h,lowerlim,upperlim,default,isint);
%
% h is an the object handle. If the string of h is not a number or
% lies outside (lowerlim upperlim) then string is set to num2str(default).
% If isint == 1 it also checks that string is an integer.

if nargin == 4
    isint = 0;
end

entry = str2double(get(h,'String'));
if ~isint
    if isnan(entry) | entry <= lowerlim | entry >= upperlim  
        errordlg(sprintf('Entry must be a number between %1.0f and %1.0f',lowerlim,upperlim),'Invalid Input','modal')
        set(h,'String',num2str(default));
    end
else
    if isnan(entry) |floor(entry)~=entry| entry < lowerlim | entry > upperlim  
        errordlg(sprintf('Entry must be an integer between %1.0f and %1.0f',lowerlim,upperlim),'Invalid Input','modal')
        set(h,'String',num2str(default));
    end
end 
    

    

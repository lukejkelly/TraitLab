function strout = strtrunc(strin,n);

% strout = strtrunc(strin,n);
% Returns the last n characters (including whitespace) of strin

if n <= 0 | n~=floor(n)
    strout = '';
    return
end
if length(strin)<=n
    strout=strin;
else
    strout=strin(end-n+1:end);
end

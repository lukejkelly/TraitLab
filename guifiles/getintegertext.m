function [x,ok] = getintegertext(h,namestr,maxx)

str = get(h,'String');
x = [];
ok = 1;
if ~isempty(str)
    x = sort(unique(str2num(str)));
    if isempty(x)
        x =[];
        ok = 0;
        errordlg([str ' is not a valid vector of ' namestr],'Invalid Input','modal')
    else
        if ~exist('maxx')
            if any(floor(x)~=x) | x(1) < 1
                x=[];
                ok = 0;
                errordlg(sprintf([namestr ' must be a vector of integers >= 1']),'Invalid Input','modal')
            end
        else
            if any(floor(x)~=x) | x(1) < 1 | x(end) > maxx
                x=[];
                ok = 0;
                errordlg(sprintf([namestr ' must be a vector of integers between 1 and %1.0f'],maxx),'Invalid Input','modal')
            end
        end
    end
end


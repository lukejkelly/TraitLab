function ok=writehtml(name, text, fig)
global IMGSIZE
% function ok=writehtml(name, text, fig) 

ok=1;

if exist([name '.html'],'file')
    %append
    fid = fopen([name '.html'],'a');
else
    % write new
    fid = fopen([name '.html'],'w');
    fprintf(fid,'<HTML><HEAD></HEAD><BODY>\n');
end

name=strrep(name,'\','\\'); % needed for Windows file names RJR 2011-03-03
text=strrep(text,'\','\\');

if nargin>=2, fprintf(fid,['<p>' text]); end

if nargin>=3
    %save figure
    figure(fig);
    suffix=[datestr(now,'yyyymmddHHMMSSFFF') int2str(fig)];
    print('-dbmp', strcat('-f',int2str(fig)), strcat(name,'fig',suffix,'.bmp'));
    
    %include figure
    fprintf(fid,['<img src=''' strcat(name,'fig',suffix,'.bmp''') ' width=''' int2str(IMGSIZE) ''' /><br/>']);
end

if nargin>=2, fprintf(fid,'</p>'); end

fprintf(fid,'\n');
fclose(fid);
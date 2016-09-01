function [content] = ObserveData(content,lang_mask,col_mask,lost)
% function [content] = ObserveData(content,lang_mask,col_mask,lost)
% masks rows lang_mask and columns col_mask; throws out data that didn't
% survive into more than lost leaves

global IN

[NS,L]=size(content.array);
if (~isempty(content.NS) & NS~=content.NS) | (~isempty(content.L) & L~=content.L)
   disp('data array passed into ObserveData with incorrect size data');keyboard;pause;
end

disp(sprintf('%g x %g data array passed into ObserveData',NS,L));

if ~isempty(col_mask)
    Ld=length(col_mask);
    disp(sprintf('Masking %d trait columns from the analysis',Ld));
    content.array(:,col_mask)=[];
end

if ~isempty(lang_mask)
    masked = char(content.language(lang_mask));
    masked = [masked char(10*ones(length(lang_mask),1))]';
    disp(sprintf('Masking the following taxa out of the analysis:%s\b',masked));
    content.array(lang_mask,:)=[];
    if ~isempty(content.language)
        content.language=CelDel(content.language,lang_mask);
    end
    disp(sprintf('Throwing out traits present in %g remaining taxa or less\n (missing data are not "present")',lost));
else
    if lost>=0
        disp(sprintf('Throwing out traits present in %g taxa or less\n (missing data are not "present")',lost));
    end
end

del=find(sum([content.array==IN])<=lost);
content.array(:,del)=[];
if ~isempty(content.cognate)
   content.cognate=CelDel(content.cognate,union(del,col_mask));
end

[content.NS,content.L]=size(content.array);

disp(sprintf('%g traits removed from data, %g remain',L-content.L,content.L));

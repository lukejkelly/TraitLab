function [ok,IC,CL,CLM,CLAM]=initclades(clademaskstr,IC,CL,cladeagesmaskstr)

if nargin <=3, cladeagesmaskstr=''; end

ok = 1;
CLM = [];
CLAM=[];

if ~isempty(cladeagesmaskstr)
    cladeagesmask = sort(unique(str2num(cladeagesmaskstr)));
    CLAM=cladeagesmask;
    if isempty(cladeagesmask)
        % str2num could not interpret vector string
        disp([cladeagesmaskstr ' is not a valid vector of clades to omit'])
        ok = 0;
    else
        % check that all numbers are natural and less than number of clades
        if any(floor(cladeagesmask)~=cladeagesmask) | cladeagesmask(1) < 1 | cladeagesmask(end) > length(CL)
            disp(sprintf('\nClades to omit must be a vector of integers between 1 and %1.0f\n',length(CL)))
            ok = 0;
        else
            keepcladeages = ones(1,length(CL));
            keepcladeages(cladeagesmask) = 0;
            disp(sprintf('\nIgnoring the age bounds for the following clades in analysis:'))
            for i = 1:length(cladeagesmask)
                disp(CL{cladeagesmask(i)}.name)
                CL{cladeagesmask(i)}.rootrange=[1 Inf];
                CL{cladeagesmask(i)}.adamrange=[];
            end
        end
    end
end
            

if ~isempty(clademaskstr)
    clademask = sort(unique(str2num(clademaskstr)));
    CLM = clademask;
    if isempty(clademask)
        % str2num could not interpret vector string
        disp([clademaskstr ' is not a valid vector of clades to omit'])
        ok = 0;
    else
        % check that all numbers are natural and less than number of clades
        if any(floor(clademask)~=clademask) | clademask(1) < 1 | clademask(end) > length(CL)
            disp(sprintf('\nClades to omit must be a vector of integers between 1 and %1.0f\n',length(CL)))
            ok = 0;
        else
            keepclade = ones(1,length(CL));
            keepclade(clademask) = 0;
            disp(sprintf('\nIgnoring the following clades in analysis:'))
            for i = 1:length(clademask)
                disp(CL{clademask(i)}.name)
            end
            if length(clademask)==length(CL)
                disp(sprintf('\nAll clades have been dropped, Clading switched off'));
                CL=[];
                IC=OFF;
            else
                disp('')
                CL = CL(logical(keepclade));
            end
        end
    end
end

if ok
    disp(sprintf('\nImposing clades:'))
    for i = 1:length(CL)
        disp(CL{i}.name)
    end
    disp('')
end
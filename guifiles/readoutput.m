function [opt,ok] = readoutput(optfile)

% [opt,ok] = readoutput(optfile) reads from the output files
% optfile.nex and optfile.txt to recover the output structure opt
% note that optfile a string with the pathname and extensionless filename
% ok == 0 if there is an error
%
% optfile.nex is a nexus file with a tree block containing Nsamp trees
% optfile.txt is a text file with 3 header lines followed by Nsamp lines each with
% space seperated fields sample number, log prior, log llkhd, root time and mu

%11/07/07 - Robin Ryder modified now 11/12 output columns was 7/8
% 22/03/14 Luke - no need for modification to account for beta as function
% can already handle a 12th column in the output file.

ok = 1;
opt = pop('output');

pro=0;
if pro
    profile on;
end

if ~exist([optfile '.nex'],'file') | ~exist([optfile '.txt'],'file')
    disp(['Error in readoutput: ' optfile '.txt or .nex not found']) 
    ok = 0;
end
    
if ok
    fid = fopen([ optfile '.txt']);
    allstats = textscan(fid,'%f','headerlines',3);
    fclose(fid);
    if isempty(allstats)
        disp(['Error in readoutput: ' optfile '.txt contains no information']) 
        ok = 0;
    end
end

if ok
    allstats = allstats{1};
    % try to determine how many columns we have
    fid = fopen([ optfile '.txt']);
    xx=textscan(fid,'%s%[^\n]','delimiter','\n','headerlines',3);
    fclose(fid); % LUKE 11/11/2015
    cols = length(str2num(xx{1}{1}));
    %expect 11 or 12 columns
    if cols < 11 || cols > 12
        disp(sprintf('Error in read output - appears there are %g columns when 11 or 12 expected.',cols))
        ok = 0;
    end
end

if ok
    Nsamp = length(allstats)/cols;
    if Nsamp ~= floor(length(allstats)/cols)
        disp(['Error in readoutput: ' optfile '.txt does not have equal length columns.'])
        ok =0;
    else
        allstats = reshape(allstats,cols,Nsamp);
    end
end

if ok
    if any(allstats(1,:)~=1:Nsamp)
        % first column is not all sample numbers
        disp([sprintf('\n') 'Error in readoutput: ' optfile '.txt does not have sample numbers as its first column.'])
        prob = min(find(allstats(1,:)~=1:Nsamp));
        disp(sprintf('Problem first found in row %1.0f\n', prob))
        ok =0;
    else  
        % file seems ok - intialise opt structure and bung in stats
        opt.stats = zeros(cols-1,Nsamp);
        if cols == 11
           opt.stats([1:5,7:11],:) = allstats(2:11,:);
        else
            opt.stats([1:5,7:12],:) = allstats(2:12,:);
        end
        opt.Nsamp = Nsamp;
    end
end

% read in the trees
if ok 
    alltrees = readalltrees([optfile '.nex']);
    ntrees = length(alltrees);
    if ntrees ~= Nsamp
        disp(sprintf('\nError in readoutput: %s.nex contains %1.0f trees while \nthe .txt file contains %1.0f samples\n',optfile,ntrees,Nsamp))
        ok =0;
    end
end

if ok
    opt.trees = alltrees;
end

% read in the cat trees
if ok && exist([optfile 'cat.nex'],'file')
    allcattrees=readalltrees([optfile 'cat.nex']);
    ncattrees=length(allcattrees);
    if ncattrees ~= Nsamp
        disp(sprintf('\nError in readoutput: %scat.nex contains %1.0f trees while \nthe .txt file contains %1.0f samples\n',optfile,ncattrees,Nsamp))
        ok=0;
    end
end

if ok && exist([optfile 'cat.nex'],'file')
    opt.cattrees=allcattrees;
    global MCMCCAT
    MCMCCAT=1;
end

if pro
    profile report readout
end
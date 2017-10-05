function llkd = converttofullLLKD(datafile,optfile,savefile)
% llkd = converttofullLLKD(datafile,optfile,savefile)
% use to calculate correct log likelihood values from old traitlab output
% give datafile with .nex extension, optfile without any extension 
% and savefile the name of file to write LLKD values to
% example: converttofullLLKD('synthdata.nex','tloutput','LLKDvalsfortloutput.txt')

global ROOT

% read in data file
[s,content] = nexus2stype(datafile);

% read in output file
opt = readoutput(optfile);

% we aren't worried about the prior but need it in makestate - assign
% arbitrary values
prior.type = 1;
prior.isclade = 0;

llkd = zeros(1,opt.Nsamp);
% step through output making tree and calcualting loglkd
for i = 1:opt.Nsamp
    % make stree
    s = rnextree(opt.trees{i});
    % make state
    state=makestate(prior,opt.stats(4,i),opt.stats(7,i),opt.stats(5,i),opt.stats(9,i),opt.stats(8,i),content,s); 
    scat=rnextree(opt.cattrees{i}); % this should give exactly the same tree; we just need to copy the number of catastrophes
    state.cat=zeros(2*state.NS);
    for node=1:(2*state.NS)
        if node.type <ROOT
            state.cat(node)=scat(scat(node.parent)).time-scat(node).time;
        end
    end
    % calculate llkd
    llkd(i) = LogLkd(state,opt.stats(7,i));
end

if exist(savefile,'file')
    ok = input(sprintf('The nominated file (%s) already exists.\nType ''y'' to overwite and proceed. ',savefile),'s');
else
    ok = 'y';
end

if strmatch(lower(ok),'y','exact')
    % write to file
    fid = fopen(savefile,'w');
    fprintf(fid,'[Log likelihood values to accompany output file %s as calculated by converttofullLLKD %s]\n',optfile,datestr(now));
    fprintf(fid,'%15.5f\n',llkd)
    fclose(fid);
else
    disp('Overwrite cancelled, file not written.')
end
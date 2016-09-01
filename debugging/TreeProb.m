function [pbar,stdp,qbar,stdq]=TreeProb(gt,outfilename,d)

%post-processing to estimate probability a given tree
%is true tree

globalswitches;
globalvalues;

[opt,ok] = readoutput(outfilename);
steps=length(d:d:opt.Nsamp);
t4=length(gt);
p=zeros(steps,t4);
q=zeros(steps,1);
for j=1:steps
    [s,errmess]=rnextree(opt.trees{j*d});
    for t=1:t4, 
        [gtt,errmess]=rnextree(gt{t}); 
         p(j,t)=equaltrees(s,gtt);
    end
    q(j)=any([s([s([s.type]==LEAF).parent]).type]==ROOT);
    %disp(sprintf('(%d) %d %d',j,steps));
end

[pbar,stdp,taup]=stats(p',round(steps/5),{},[]);
for k=1:t4
    disp(sprintf('tree %d, prob=%5.3f(%5.3f)',k,pbar(k),stdp(k)));
end

[qbar,stdq,tauq]=stats(q',round(steps/5),{},[]);
disp(sprintf('Any tree with single-leaf outgroup, prob=%5.3f(%5.3f)',qbar,stdq));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check uniform topo (and Markov move 2) on 4 leaf tree
%check that, if topo is off, get weighting by number of linear extensions
%There was also a problem with move 2 NNI in Markov (Echoose.m) - it had
%a while loop that kept trying for a legal swap - that was wrong
%without a correcting logqq - GKN 1/4/11

%to check each of the topo-changing mcmc updates individualy set 
%move=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%in fullsetup.m (etc)  

%Either - load synthdata.nex and run it with traits 1:29 removed 
%(and cats/missing/vary mu) switched off 
%Or use the batch below - all the files are in this dir

%%
%might want to set verbose to 1 in batchTraitlab;
%you need really long runs for debugging - presently set short
batchTraitlab;

%%
gt={'(4,((1,2),3))','(4,((1,3),2))','(2,((1,3),4))','(3,((1,2),4))','(3,((1,4),2))','(2,((1,4),3))','(1,((2,3),4))','(4,((2,3),1))',...
     '(4,((2,3),1))', '(1,((2,3),4))','(2,((3,4),1))','(1,((3,4),2))', '((1,2),(3,4))','((1,3),(2,4))','((1,4),(2,3))'};

%short - qr is the probability to get an unbalanced tree
%qr=2/3 if UT off, qr=12/15 if UT on
%[pr,er,qr,eq]=TreeProb({gt{[1,15]}},'tloutput',1);

%long
[pr,er,qr,eq]=TreeProb(gt,'tloutput',1);
%should see [12 1/18's and 3 1/9's] if uniform topo is off
%should see 15 1/15's if UT is on
figure; errorbar(pr,2*er,'x')
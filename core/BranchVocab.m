function [nextwords,nextL]=BranchVocab(prevwords,p)

prevL=length(prevwords);
nextwords=prevwords(rand(1,prevL)<p);
nextL=length(nextwords);
function synthetizeinblocks(s)
% synthetize data in blocks
% no catastrophes
% assumes all leaves are modern
% assumes
% RJR 2009-08-09

global LEAF ROOT

leaves=find([s.type]==LEAF);
NS=length(leaves);
Root=find([s.type]==ROOT);

nmc=100; %number of meaning categories
tpc=1.3; %average number of traits per category in a given language

tnt=0; %total number of traits

ListOfTimes=[s.time];
ListOfTimes=sort(ListOfTimes(ListOfTimes>0),1,'descend'); % take out leaves


for i=1:nmc
     
    nw=MyPoisson(tpc);
    data{
    for j=2:NS
        
function array=SimMis(array,s,L,L_vals)
% some data go missing
% RJR 29/05/08
% if L_vals is defined then data will go missing in blocks

global LEAF MIX

leaves=find([s.type]==LEAF);
NS=length(leaves);

if nargin<4, L_vals=ones(1,L); end
nmc=length(L_vals);
Lcum=[0,cumsum(L_vals)];

for k=1:NS
    missing=find(rand(nmc,1)>s(leaves(k)).xi); %was s(k).xi but ??? GKN 1/4/11
    nmiss=length(missing);
    if nmiss==nmc
        disp(['Warning: all data for leaf ' num2str(k) ' are missing.']);
    end
    for j=1:nmiss
        array(k,(Lcum(missing(j))+1):Lcum(missing(j)+1))=MIX;
    end
end


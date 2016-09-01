 function s = insp(output,data,pos,cog)

GlobalSwitches;
global za zb zc;
% merge data and tree
s = rnextree(output.trees{pos});
[s,ok]=mergetreedata(s,data.array,data.language);
if ok     
    % merge work, now need to fill work variables in tree
    state = pop('state');
    state.NS = length(s)/2+1;
    state.L = size(data.array,2);
    state.tree =s;
    state.root = find([s.type] == ROOT);
    state.leaves = find([s.type] == LEAF);
    state.nodes = find([s.type]==ANST);
    [state.tree([state.leaves state.nodes]).mark]=deal(1);
    za=zeros(1,state.L);
    zb=zeros(1,state.L);
    zc=zeros(1,state.L);
    [state.tree]=ActiveI(state.tree,state.root);
    state.tree(state.root).CovI=[state.tree(state.root).ActI{:}];
    state.tree=CoversI(state.tree,state.root);
    % get drawing options
    %vals = get([handles.leafnamesrb handles.anstnumsrb handles.showcovrb handles.showcogcb],'Value');    
    vals={0,1,0,0};
    if vals{1}
        vb = LEAF;
    elseif vals{2}
        vb = ANST;
    elseif vals{3}
        vb = COGS;
    else
        vb = -1;
    end
    showcog = 0;
    cogname = '';
    if vals{4}
        showcog = cog;
        if length(handles.data.cognate)>=showcog
            cogname = data.cognate{showcog};
        end
    end
    draw(state.tree,output.treefig,vb,['Sample number ' num2str(pos)],showcog,cogname);
    s= state.tree;
end

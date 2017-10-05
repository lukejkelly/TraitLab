function output=getcladetimes(data,output,clade)

nt=size(output.trees,1);
if output.Nsamp~=nt, disp('output.Nsamp~=size(output.trees,1)'); keyboard;pause; end

nc=size(clade,2);

mu=output.stats(4,:);
p=output.stats(6,:);

cladetimes=zeros(nc,nt);

for k=1:nt
    s=rnextree(output.trees{k});
    state=tree2state(s);
    state.tree = treeclades(state.tree,clade);   
    state=UpdateClades(state,[state.leaves,state.nodes],nc);
    ci=find(state.claderoot~=0);
    cladetimes(ci,k)=[state.tree([state.claderoot(ci)]).time]';
end

output.cladetimes=cladetimes;

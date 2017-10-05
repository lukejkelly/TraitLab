function result=check(state,stateTR)

%global ROOT ADAM LEAF ANST DONTMOVECATS
GlobalSwitches;
global BORROWING

[NS,N,L,mu,s,Root,LEAVES,NODES,claderoot,cllkd,olp]=new2old(state);

%if ~isempty(stateTR)
%   [NStr,Ntr,Ltr,muTR,sTR,RootTR,LEAVEStr,NODEStr,claderootTR,llkdTR,lpTR]=new2old(stateTR);
%end

result=[];

for k=1:N
  if s(k).type==ROOT
     if k~=Root
       result=[result,[1.1;k]];
     end
     if s(s(k).parent).type~=ADAM
       result=[result,[1.2;k]];
     end
     if s(s(k).parent).child~=Root
       result=[result,[1.3;k]];
     end
     %if s(s(k).parent).time~=realmax
     %   result=[result,[1.4;k]];
     %end
  elseif s(k).type==LEAF
     if ~isempty(s(k).child)
        result=[result,[2.2;k]];
     end
     if s(s(k).parent).child(s(k).sibling)~=k
        result=[result,[2;k]];
     end
     if s(s(k).parent).time<s(k).time
        result=[result,[2.1;k]];
     end
  else
     if s(k).type~=ANST
         result=[result,[3.1;k]];
     end
     if s(s(k).parent).child(s(k).sibling)~=k
        result=[result,[3.2;k]];
     end
     if s(s(k).parent).time<s(k).time
        result=[result,[3.3;k]];
     end
  end
end

for j=[LEAVES,NODES]
   if s(j).mark~=0
		 result=[result,[4;j]];
   end
end

%s=WorkVars(NS,L,s,Root); %unused assignment RJR 05/07/07
TOPOLOGY=1;
state=MarkRcurs(state,[LEAVES,NODES],TOPOLOGY); 
%TODO compare c and state.claderoot

%llkd=LogLkd(state); % Luke --- moved this inside if/else.
if BORROWING 
    [llkd, ~] = logLkd2(state);
else
    llkd=LogLkd(state);
end

if (L>0) && (abs(cllkd-llkd)/abs(llkd) > 1e-10) % Luke --- this may fail when doing EA method.
    cllkd
    llkd
   result=[result,[5;Root]];
end

% Check the number of catastrophes
if state.ncat ~= sum(state.cat)
    result=[result, [6;sum(state.cat)]];
end


if DONTMOVECATS
    for k=1:N
        if s(k).type<=ANST && state.cat(k)~=1, result=[result, [7;k]]; end
    end
end
    

% Check that the length is correct
if (state.length-TreeLength(state.tree,state.root))/state.length>1e-10
    result=[result, [8;state.length]];
end

if isempty(result)
  result=0;
end

function [i,newage,logq]=Schoose(state)
%function [i,newage,logq]=Schoose(state)

global ROOT

s=state.tree;
%if rand<0.5
%   t=s(state.nodes).time;
%   v=cumsum(t);
%   tot=v(2*state.NS-1);
%   r=rand*tot;
%   i=state.nodes(min(find(r<v)));
%else
   i=state.nodes(ceil(rand*(state.NS-1)));
   %end

iT=s(i).time;
k=s(i).parent; 
j1=s(i).child(1);
j2=s(i).child(2);
jT=max(s(j1).time,s(j2).time);

if s(i).type==ROOT
  %this way to update root should be more adaptive to likelihood
  tau=iT-jT;
  delta=(-0.5+rand*1.5)*tau;
  taup=tau+delta;
  logq=log(tau/taup);
  newage=jT+taup;
else
  newage=jT+rand*(s(k).time-jT);
  logq=0;
end

function [nextstate,Accept]=MHupd(model,state,nstate,logq,update);

GlobalSwitches;

logr=nstate.logprior-state.logprior+nstate.loglkd-state.loglkd+logq;

if ( (logr>0) | (log(rand)<logr) )
   nextstate=nstate;
   Accept=1;
else
   nextstate=state;
   Accept=0;
end

if TESTUP & check(state,[])
   disp(['Error in update:',update]);
   keyboard;pause;
end
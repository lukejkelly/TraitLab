function seq=GTRsim(s,k,seq,L,mu,Q)

for kc=s(k).child
  t=(s(k).time-s(kc).time)*mu; % *mu not /mu since "time" is yrs and we want mutns
  eQt=expm(Q*t);
  for site=1:L
    seq(kc,site)=disample(eQt(seq(k,site),:));
  end
  %keyboard;
  seq=GTRsim(s,kc,seq,L,mu,Q);
end


function B=IsClade(v,s,r)

a=mrca(v,s,r);
m=NumLeavesUnder(a,s);
if (m<length(v)), error('IsClade found a MRCA that doesnt cover its leaves'); end
B=(m==length(v));
function s=WorkVars(NS,L,s,Root)
%function s=WorkVars(NS,L,s,Root)

[s.ActI]=deal({[]  []  []});
[s.CovI]=deal([]);
[s.difCovI]=deal([]);
%[s.Tw]=deal(zeros(1,NS));
[s.w]=deal(zeros(1,L));
[s.PD]=deal(zeros(1,L));
[s.PDcat]=deal(zeros(1,L));

s(Root).CovI=1:L;
s(s(Root).parent).CovI=1:L;

s(Root).TCovI=1:NS;
s(s(Root).parent).TCovI=1:NS;

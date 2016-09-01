function clade=claderootrangetimes(clade,i)

%function clade=claderootrangetimes(clade,topclade)

for m=clade{i}.children
    clade{m}.length=clade{i}.rootrange(1)-clade{m}.rootrange(1);
    clade=claderootrangetimes(clade,m);
end
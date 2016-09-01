function lcaz=CelDel(caz,langoutz)

nsz=size(caz,1);
langsz=1:nsz;
langinz=setdiff(langsz,langoutz);
lnsz=length(langinz);
lcaz=cell(lnsz,1);
if ~isempty(lcaz)
   [lcaz{:}]=deal(caz{langinz});
end
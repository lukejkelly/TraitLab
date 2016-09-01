function t=ifc(s,NODES)

[dump1,ifc]=min([s(NODES).time]);
t=s(NODES(ifc)).time;

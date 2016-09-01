function [HPDlow, HPDhigh]=hpd(v,a)
% function [min max]=hpd(v,a)
% Computes a (100-a)% HPD for the vector v
% assumes a unimodal distribution
% Robin J. Ryder, 02/10/2008

v=sort(v);
l=length(v);
n=round(l*a/100)+1;

intL=v(l-n+1:l)-v(1:n);
[foo, i]=min(intL);

HPDlow=v(i);
HPDhigh=v(l-n+i);
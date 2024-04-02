function d=tdiss(x,y)
%dissipation function

d=(x-y).*log(x./y);
d(x<=0)=0;
d(y<=0)=0;


function F=Goodwin(x,k,d,n)
% the original Goodwin model itself
F=zeros(3,1);
F(1)=k(1)/(1+x(3)^n)-d(1)*x(1);
F(2)=k(2)*x(1)-d(2)*x(2);
F(3)=k(3)*x(2)-d(3)*x(3);

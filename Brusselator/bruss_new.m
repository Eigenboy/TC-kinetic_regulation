function f=bruss_new(x,k,a)
%brusselator itself
f=zeros(2,1);
f(1)=(k(1,1)+k(1,8))*a-(k(1,5)+k(1,4)+k(1,2))*x(1)+k(1,6)*x(2)+k(1,3)*x(1)^2*x(2)-k(1,7)*x(1)^3;
f(2)=k(1,2)*x(1)-k(1,6)*x(2)-k(1,3)*x(1)^2*x(2)+k(1,7)*x(1)^3;

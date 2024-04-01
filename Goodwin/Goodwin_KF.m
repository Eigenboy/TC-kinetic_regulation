function F=Goodwin_KF(x,k,d,A,Kd)
% the modified Goodwin model proposed by Kim and Forger
F=zeros(3,1);
fun = (A - x(3) - Kd + sqrt((A - x(3) - Kd)^2 + 4*A*Kd))/(2*A);
F(1)=k(1)*fun-d(1)*x(1);
F(2)=k(2)*x(1)-d(2)*x(2);
F(3)=k(3)*x(2)-d(3)*x(3);

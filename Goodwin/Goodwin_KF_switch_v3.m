function F=Goodwin_KF_switch_v3(x,k,d,A,Kd,Km,Km2)
% the simplified Zhou-Kim model itself
F=zeros(4,1);
fun = (A - x(4) - Kd + sqrt((A - x(4) - Kd)^2 + 4*A*Kd))/(2*A);
F(1) = k(1)*fun-d(1)*x(1);
F(2) = k(2)*x(1) - k(3)*x(2)/(x(2)+Km) - k(4)*x(2)/(x(2)+Km2) - d(2)*x(2);
F(3) = k(3)*x(2)/(x(2)+Km) - d(3)*x(3);
F(4) = k(4)*x(2)/(x(2)+Km2) - d(4)*x(4);
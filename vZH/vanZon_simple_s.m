function dC=vanZon_simple_s(C,kp,kdp,kf,kb,f6,b0,At,Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s)
% The van Zon - Hatakeyama model iself

dC=zeros(20,1);
C21=max(0,Ct-sum(C));%inactive C6
A=max(0,At-(C(8,1)+C(9,1)+C(10,1)+C(11,1)+C(12,1)+C(13,1)+C(14,1)));

%active KaiC 1~7
dC(1)=b0*C(15,1)-(gamma_tb*b0+kf(1)*A)*C(1,1)+kb(1)*C(8,1);
dC(2)=-kf(2)*C(2,1)*A+kb(2)*C(9,1);
dC(3)=-kf(3)*C(3,1)*A+kb(3)*C(10,1);
dC(4)=-kf(4)*C(4,1)*A+kb(4)*C(11,1);
dC(5)=-kf(5)*C(5,1)*A+kb(5)*C(12,1);
dC(6)=-kf(6)*C(6,1)*A+kb(6)*C(13,1);
dC(7)=gamma_tf*f6*C21-(f6+kf(7)*A)*C(7,1)+kb(7)*C(14,1);

%KaiAKaiC 1~7
dC(8) =kf(1)*A*C(1,1)-(kb(1)+alpha(1)*kp(1))*C(8,1)+gamma*alpha(1)*kp(1)*C(9,1);
dC(9) =kf(2)*A*C(2,1)-(kb(2)+alpha(2)*kp(2)+gamma*alpha(1)*kp(1))*C(9,1)+alpha(1)*kp(1)*C(8,1)+gamma*alpha(2)*kp(2)*C(10,1);
dC(10)=kf(3)*A*C(3,1)-(kb(3)+alpha(3)*kp(3)+gamma*alpha(2)*kp(2))*C(10,1)+alpha(2)*kp(2)*C(9,1)+gamma*alpha(3)*kp(3)*C(11,1);
dC(11)=kf(4)*A*C(4,1)-(kb(4)+alpha(4)*kp(4)+gamma*alpha(3)*kp(3))*C(11,1)+alpha(3)*kp(3)*C(10,1)+gamma*alpha(4)*kp(4)*C(12,1);
dC(12)=kf(5)*A*C(5,1)-(kb(5)+alpha(5)*kp(5)+gamma*alpha(4)*kp(4))*C(12,1)+alpha(4)*kp(4)*C(11,1)+gamma*alpha(5)*kp(5)*C(13,1);
dC(13)=kf(6)*A*C(6,1)-(kb(6)+alpha(6)*kp(6)+gamma*alpha(5)*kp(5))*C(13,1)+alpha(5)*kp(5)*C(12,1)+gamma_s*alpha(6)*kp(6)*C(14,1);
dC(14)=kf(7)*A*C(7,1)-(kb(7)+gamma_s*alpha(6)*kp(6))*C(14,1)+alpha(6)*kp(6)*C(13,1);

%inactive KaiC 1~7
dC(15)=gamma_tb*b0*C(1,1)+kdp(1)*C(16,1)-(b0+gamma2*kdp(1))*C(15,1);
dC(16)=gamma2*kdp(1)*C(15,1)+kdp(2)*C(17,1)-(kdp(1)+gamma2*kdp(2))*C(16,1);
dC(17)=gamma2*kdp(2)*C(16,1)+kdp(3)*C(18,1)-(kdp(2)+gamma2*kdp(3))*C(17,1);
dC(18)=gamma2*kdp(3)*C(17,1)+kdp(4)*C(19,1)-((kdp(3)+gamma2*kdp(4)))*C(18,1);
dC(19)=gamma2*kdp(4)*C(18,1)+kdp(5)*C(20,1)-((kdp(4)+gamma2*kdp(5)))*C(19,1);
dC(20)=gamma2*kdp(5)*C(19,1)+kdp(6)*C21-(kdp(5)+gamma2*kdp(6))*C(20,1);

% dC(21)=gamma2*kdp(6)*C(20,1)+f6*C(7,1)-(kdp(6)+gamma_tf*f6)*C21;
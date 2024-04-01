function f=Tyson_per_2D(Z,v_m,k_m,v_p,k_1,k_2,k_3,J,K,A)
% Tyson's 2D model for Drosophila circadian clock

f = zeros(2,1);
q = 2/(1+sqrt(1+8*K*Z(2)));
P1 = q*Z(2);
P2 = (1-q)*Z(2)/2;
f(1) = v_m/(1+P2^2/A^2) - k_m*Z(1);
f(2) = v_p*Z(1) -(k_1*P1+2*k_2*P2)/(J+Z(2)) - k_3*Z(2);
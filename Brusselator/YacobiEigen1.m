function lambda=YacobiEigen1(k,x_s,y_s)
%Eigenvalues of the Jacobian matrix of Brusselator
J=[-(k(1,4)+k(1,2)+k(1,5))+2*k(1,3)*x_s*y_s-3*k(1,7)*x_s^2,k(1,6)+k(1,3)*x_s^2;
    k(1,2)-2*k(1,3)*x_s*y_s+3*k(1,7)*x_s^2,-(k(1,6)+k(1,3)*x_s^2)];
lambda=eig(J);

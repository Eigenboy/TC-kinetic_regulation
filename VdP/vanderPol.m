function F=vanderPol(x,mu,w)
% the van der Pol (VdP) model itself

F=zeros(2,1);
F(1)=-w*x(2)-mu*x(1)*(x(1)^2/3-1);
F(2)=w*x(1);
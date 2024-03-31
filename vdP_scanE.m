% VdP model; fixed μ (mu0) and Ω (w);
% scan the activation-energy space from 10 kT to 30 kT;

clear;
tic
w0=1;
mu0=30;
T=[298;308];

[t1,X1]=ode45(@(t,x) vanderPol(x,mu0,w0),[0,50*pi],[0.01,0.01]);
[t1,X1]=ode45(@(t,x) vanderPol(x,mu0,w0),[0,250*pi],X1(end,:));
P1=periods(t1,X1(:,1));
%plot(X1(:,1),X1(:,2));
N=101;
hwait=waitbar(0,'Waiting');
for ii=1:N
    for jj=1:N
        waitbar(ii/N,hwait,'Waiting');
        E1(ii,jj)=(9.8+0.2*ii)*8.31*T(1);
        E2(ii,jj)=(9.8+0.2*jj)*8.31*T(1);

        mu=mu0*exp(E1(ii,jj)/(8.31*T(1)))*exp(-E1(ii,jj)/(8.31*T(2)));
        w=w0*exp(E2(ii,jj)/(8.31*T(1)))*exp(-E2(ii,jj)/(8.31*T(2)));
    
        [t2,X2]=ode45(@(t,x) vanderPol(x,mu,w),[0,50*pi],[0.01,0.01]);
        [t2,X2]=ode45(@(t,x) vanderPol(x,mu,w),[0,250*pi],X2(end,:));

        P2=periods(t2,X2(:,1));
        Q10(ii,jj)=P2/P1;
    
    end
end
toc

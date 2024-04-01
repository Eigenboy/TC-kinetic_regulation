% VdP model; vary μ (mu0) and keep Ω (w) unchanged;
% scan the activation-energy space and compute the TC domain (%) by q/N^2;

clear;
w0=1;
T=[298;308];

M = 41;
for kk=1:M
    for jj=1:M
        E1(kk,jj)=(9.5+0.5*kk)*8.31*T(1);
        E2(kk,jj)=(9.5+0.5*jj)*8.31*T(1);
    end
end

N = 31;
hwait=waitbar(0,'Waiting');
for ii=1:N
    waitbar(ii/N,hwait,'Waiting');
    mu0(ii)=1*ii-1;

    [t1,X1]=ode45(@(t,x) vanderPol(x,mu0(ii),w0),[0,50*pi],[0.01,0.01]);
    [t1,X1]=ode45(@(t,x) vanderPol(x,mu0(ii),w0),[0,250*pi],X1(end,:));
    %plot(X1(:,1),X1(:,2));
    P1(ii)=periods(t1,X1(:,1));
    q(ii)=0;
    for kk=1:M
        for jj=1:M
            mu=mu0(ii)*exp(E1(kk,jj)/(8.31*T(1)))*exp(-E1(kk,jj)/(8.31*T(2)));
            w=w0*exp(E2(kk,jj)/(8.31*T(1)))*exp(-E2(kk,jj)/(8.31*T(2)));

            [t2,X2]=ode45(@(t,x) vanderPol(x,mu,w),[0,50*pi],[0.01,0.01]);
            [t2,X2]=ode45(@(t,x) vanderPol(x,mu,w),[0,250*pi],X2(end,:));
            P2=periods(t2,X2(:,1));
            Q10=P2/P1(ii);
            if Q10>0.9 && Q10<1.1
                q(ii)=q(ii)+1;
            end
            
        end
    end
end

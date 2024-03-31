% VdP model; vary both μ (mu0) and Ω (w)
% fixed activation energies

clear;
tic;

T=[298;308];
E=8.31*T(1)*[13.3,26.7];

N=41;
M=41;
hwait=waitbar(0,'Waiting');
for ii=1:N
    waitbar(ii/N,hwait,'Waiting');
    for jj = 1:M
        mu0(ii,jj) = 10^(-1.1 + 0.1*ii);
        w0(ii,jj) = 10^(-1.1 + 0.1*jj);

        w=[w0(ii,jj),w0(ii,jj)*exp(E(1)/(8.31*T(1)))*exp(-E(1)/(8.31*T(2)))];
        mu=[mu0(ii,jj),mu0(ii,jj)*exp(E(2)/(8.31*T(1)))*exp(-E(2)/(8.31*T(2)))];
        P_estimate = max((3-2*log(2))*mu0(ii,jj)/w0(ii,jj)^2,2*pi/w0(ii,jj));
    
        [t1,X1]=ode23s(@(t,x) vanderPol(x,mu(1),w(1)),[0,25*P_estimate],[0.01,0.01]);
        [t1,X1]=ode23s(@(t,x) vanderPol(x,mu(1),w(1)),[0,50*P_estimate],X1(end,:));
        %plot(X1(:,1),X1(:,2));
        
        [t2,X2]=ode23s(@(t,x) vanderPol(x,mu(2),w(2)),[0,25*P_estimate],[0.01,0.01]);
        [t2,X2]=ode23s(@(t,x) vanderPol(x,mu(2),w(2)),[0,50*P_estimate],X2(end,:));
    
        P1(ii,jj)=periods(t1,X1(:,1));
        P2(ii,jj)=periods(t2,X2(:,1));
        Q10(ii,jj)=P2(ii,jj)/P1(ii,jj);
        
    end
end
toc;
% VdP model; vary μ (mu0) and keep Ω (w) unchanged
% fixed activation energies

clear;

% Ω value at T=298K
w0 = 1;

% two temperatures
T = [298; 308];

% activation energies
E = 8.31 * T(1) * [13.3, 26.7];

% Ω values at two temperatures
w = [w0, w0 * exp(E(1)/(8.31*T(1))) * exp(-E(1)/(8.31*T(2)))];

N = 41;
hwait = waitbar(0,'Waiting');
for ii=1:N
    waitbar(ii/N,hwait,'Waiting');

    % μ value at T=298K
    mu0(ii) = 10^(-1.1+0.1*ii);

    % μ value at two temperatures
    mu=[mu0(ii), mu0(ii) * exp(E(2)/(8.31*T(1))) * exp(-E(2)/(8.31*T(2)))];

    % estimation of the period
    P_estimate = max((3-2*log(2))*mu0(ii)/w0^2, 2*pi/w0);
    
    % trajectories at T=298K
    [t1,X1] = ode23s(@(t,x) vanderPol(x,mu(1),w(1)),[0,25*P_estimate],[0.01,0.01]);
    [t1,X1] = ode23s(@(t,x) vanderPol(x,mu(1),w(1)),[0,50*P_estimate],X1(end,:));
    %plot(X1(:,1),X1(:,2));
    
    % trajectories at T=308K
    [t2,X2] = ode23s(@(t,x) vanderPol(x,mu(2),w(2)),[0,25*P_estimate],[0.01,0.01]);
    [t2,X2] = ode23s(@(t,x) vanderPol(x,mu(2),w(2)),[0,50*P_estimate],X2(end,:));
    
    % computing periods and Q10
    P1(ii) = periods(t1,X1(:,1));
    P2(ii) = periods(t2,X2(:,1));
    Q10(ii) = P2(ii)/P1(ii);
    
    % computing period sensitivities C_mu and C_w
    amp = 0.01;
    mu_vary = mu0(ii)*(1+amp);
    [t3,X3] = ode23s(@(t,x) vanderPol(x,mu_vary,w(1)),[0,25*P_estimate],[0.01,0.01]);
    [t3,X3] = ode23s(@(t,x) vanderPol(x,mu_vary,w(1)),[0,50*P_estimate],X3(end,:));
    P3 = periods(t3,X3(:,1));
    C_mu(ii) = 1/amp*(P3/P1(ii)-1);
    
    w_vary=w(1)*(1+amp);
    [t4,X4]=ode23s(@(t,x) vanderPol(x,mu(1),w_vary),[0,25*P_estimate],[0.01,0.01]);
    [t4,X4]=ode23s(@(t,x) vanderPol(x,mu(1),w_vary),[0,50*P_estimate],X4(end,:));
    P4=periods(t4,X4(:,1));
    C_w(ii)=1/amp*(P4/P1(ii)-1);
    
end

% computing progression speeds
% for ii=1:length(X1)
%     v_x(ii,1)=-w0*X1(ii,2)-mu(1)*X1(ii,1)*(X1(ii,1)^2/3-1);
%     v_y(ii,1)=w0*X1(ii,1);
% end
% v_abs=sqrt(v_x.^2+v_y.^2);

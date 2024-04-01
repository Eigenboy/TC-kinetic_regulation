% Simplifed Zhou-Kim model; Checking oscillation phase separation by
% computing amplitude/speed/time sensitivities

clear;
tic;
T=[298;308];
A = 0.1;
Kd = 10^(-6);
Km = 1;
Km2 = 1;
Ed=[10,10,10,10]*8.31*T(1);Ek=[25,25,5,15]*8.31*T(1);
d(1,1:4)= [1,0.01,10,1];
d(2,1:4)=d(1,1:4).*exp(Ed/(8.31*T(1))).*exp(-Ed/(8.31*T(2)));

k(1,1:4)= [10,2,10,10];
k(2,1:4)=k(1,1:4).*exp(Ek/(8.31*T(1))).*exp(-Ek/(8.31*T(2)));

M = 11;
[P1,k_p,amp_fast,tau_fast,v_fast,amp_slow,tau_slow,v_slow] = deal(zeros(1,M));

hwait=waitbar(0,'waiting');
for jj = 1:M
    waitbar(jj/M,hwait,'waiting');
    k_p(1,jj) = k(1,2)*(1 + 0.01*(jj-1));
    kv(1,1:4)=k(1,1:4);
    kv(1,5:8)=d(1,1:4);
    kv(1,2) = k_p(1,jj);

    [~,X1]=ode89(@(t,x) Goodwin_KF_switch_v3(x,kv(1,1:4),kv(1,5:8),A,Kd,Km,Km2),[0,400],[1,0,0,0]);
    [t1,X1]=ode89(@(t,x) Goodwin_KF_switch_v3(x,kv(1,1:4),kv(1,5:8),A,Kd,Km,Km2),[0,400],X1(end,:));
    % figure(1);plot(t1,X1);
    % figure(2);plot(X1(:,1),X1(:,4));hold on;
    P1(1,jj) = periods(t1(floor(length(t1)/2):end,1),X1(floor(length(t1)/2):end,4));
    
    dt = 0.0005;
    tt = 0:dt:t1(end-1);
    X_fit = spline(t1,X1(:,1)+X1(:,2),tt);
    Y_fit = spline(t1,X1(:,3)+X1(:,4),tt);
    
    [pks,ind_pks] = findpeaks(X_fit);
    X_fit_last = X_fit(ind_pks(end-1):ind_pks(end));
    Y_fit_last = Y_fit(ind_pks(end-1):ind_pks(end));
    
    [trgh,ind_trgh] = findpeaks(-X_fit_last);
    trgh = - trgh;
    
    X_slow = X_fit_last(1:ind_trgh);
    Y_slow = Y_fit_last(1:ind_trgh);
    X_fast = X_fit_last(ind_trgh:end);
    Y_fast = Y_fit_last(ind_trgh:end);
    
    figure(2);plot(X_slow,Y_slow,X_fast,Y_fast);hold on;
    
    dX_fast = diff(X_fast);
    dY_fast = diff(Y_fast);
    amp_fast(jj) = sum(sqrt(dX_fast.^2+dY_fast.^2));
    tau_fast(jj) = dt*(length(X_fit_last) - ind_trgh);
    v_fast(jj) = amp_fast(jj)/tau_fast(jj);
    dX_slow = diff(X_slow);
    dY_slow = diff(Y_slow);
    amp_slow(jj) = sum(sqrt(dX_slow.^2+dY_slow.^2));
    tau_slow(jj) = dt*(ind_trgh - 1);
    v_slow(jj) = amp_slow(jj)/tau_slow(jj);
end

k_p_norm = k_p/k_p(ceil(end/2));
amp_fast_norm = amp_fast/amp_fast(ceil(end/2));
tau_fast_norm = tau_fast/tau_fast(ceil(end/2));
v_fast_norm = v_fast/v_fast(ceil(end/2));

amp_slow_norm = amp_slow/amp_slow(ceil(end/2));
tau_slow_norm = tau_slow/tau_slow(ceil(end/2));
v_slow_norm = v_slow/v_slow(ceil(end/2));

toc;

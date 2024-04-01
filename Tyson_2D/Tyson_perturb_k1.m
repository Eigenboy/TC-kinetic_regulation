% Tyson-2D model; checking oscillation phase separation

clear;

v_m = 1*1.0;
k_m = 0.1*1.0;
v_p = 0.5*1.0;
% k_1 = 10*1.0;
k_1 = exp(2.1);
k_2 = 0.03*1.0;
k_3 = 0.1*1.0;
J = 0.05;
K = 200*1;
A = 0.1;

P_temp = 30;
N = 13;
hwait=waitbar(0,'Waiting');
for kk=1:N
    waitbar(kk/N,hwait,'Waiting');
    k_1_vary(1,kk)= k_1 * (1 - 0.14 + 0.02*kk);

    [t1,z1]=ode23s(@(t,z) Tyson_per_2D(z,v_m,k_m,v_p,k_1_vary(1,kk),k_2,k_3,J,K,A),[0,30*P_temp],[1,1]);
    P_temp1 = periods(t1,z1(:,1));
    if P_temp1 > 0
        [t1,z1]=ode23s(@(t,z) Tyson_per_2D(z,v_m,k_m,v_p,k_1_vary(1,kk),k_2,k_3,J,K,A),[0,60*P_temp1],z1(end,:));
    end
    P1 = periods(t1,z1(:,1));
    % figure(1);
    % plot(z1(:,1),z1(:,2));hold on;
    
    % P_T nullcline
    pt = 0:0.01:10;
    q = 2./(1+sqrt(1+8*K*pt));
    p_1 = q.*pt;
    p_2 = (1-q).*pt/2;
    m = (k_3*pt + (k_1_vary*p_1+2*k_2*p_2)./(J+pt))/v_p;

    figure(1);
    plot(m,pt);hold on;

    % progression speed
    speed = zeros(size(z1));
    for ii = 1:length(t1)
        speed(ii,:) = Tyson_per_2D(z1(ii,:),v_m,k_m,v_p,k_1_vary,k_2,k_3,J,K,A);
    end
    speed_abs = sqrt(speed(:,1).^2+speed(:,2).^2);
    % patch([z1(end-1000:end,1)' NaN],[z1(end-1000:end,2)' NaN],[1./speed_abs(end-1000:end,1)' NaN],'EdgeColor','interp','MarkerFaceColor','flat')
    
    % separated oscillation phases
    dt = 0.05;
    tt = 0:dt:t1(end-1);
    z1_1_fit = spline(t1,z1(:,1),tt);
    z1_2_fit = spline(t1,z1(:,2),tt);

    [pks,ind_pks] = findpeaks(z1_1_fit);
    z1_1_thresh = z1_1_fit(ind_pks(end));
    [pks2,ind_pks2] = findpeaks(-z1_2_fit(ind_pks(end-1):ind_pks(end)));
    z1_2_thresh = z1_2_fit(ind_pks2);
    ind_pks_mid = ind_pks(end-1) + ind_pks2 - 1;

    figure(1);
    plot(z1_1_fit(ind_pks(end-1):ind_pks_mid),z1_2_fit(ind_pks(end-1):ind_pks_mid));
    hold on;
    plot(z1_1_fit(ind_pks_mid:ind_pks(end)),z1_2_fit(ind_pks_mid:ind_pks(end)));

    dX_fast = diff(z1_1_fit(ind_pks(end-1):ind_pks_mid));
    dY_fast = diff(z1_2_fit(ind_pks(end-1):ind_pks_mid));
    amp_fast(kk) = sum(sqrt(dX_fast.^2+dY_fast.^2));
    tau_fast(kk) = tt(ind_pks_mid) - tt(ind_pks(end-1));
    v_fast(kk) = amp_fast(kk)/tau_fast(kk);

    dX_slow = diff(z1_1_fit(ind_pks_mid:ind_pks(end)));
    dY_slow = diff(z1_2_fit(ind_pks_mid:ind_pks(end)));
    amp_slow(kk) = sum(sqrt(dX_slow.^2+dY_slow.^2));
    tau_slow(kk) = tt(ind_pks(end)) - tt(ind_pks_mid);
    v_slow(kk) = amp_slow(kk)/tau_slow(kk);
end
k_1_norm = k_1_vary/k_1;
amp_fast_norm = amp_fast/amp_fast((N+1)/2);
tau_fast_norm = tau_fast/tau_fast((N+1)/2);
v_fast_norm = v_fast/v_fast((N+1)/2);

amp_slow_norm = amp_slow/amp_slow((N+1)/2);
tau_slow_norm = tau_slow/tau_slow((N+1)/2);
v_slow_norm = v_slow/v_slow((N+1)/2);
figure;
plot(k_1_norm,amp_fast_norm,k_1_norm,amp_slow_norm);
figure
plot(k_1_norm,tau_fast_norm,k_1_norm,tau_slow_norm);
figure
plot(k_1_norm,v_fast_norm,k_1_norm,v_slow_norm);
% save('data_perturb_k1.mat','k_1_norm','amp_fast_norm','tau_fast_norm','v_fast_norm','amp_slow_norm','tau_slow_norm','v_slow_norm','v_m','k_m','v_p','k_2','k_3','J','K','A');
% van Zon - Hatakeyama model; exploring oscillation phase separation at
% certain total KaiA

clear;
tic
fctr=0.333;
kp_0(1:6)=3.1*fctr;
kp_0(5)=1.4*kp_0(5);
kp_0(6)=1.2*kp_0(6);
kdp_0(1:6)=9.5*fctr;
b0_0=1*fctr;
f6_0=1*fctr;
kf_0(1:7)=5*1e7*fctr;
kb_0=zeros(1,6);
strgth=10;
for jj=1:7
    kb_0(jj)=strgth^7*strgth^(jj-7)*fctr;
end
Ct = 11.4;
At = 1.5;
% At = 4.2;

gamma0=exp(-8);
gamma1=exp(-8);%ADP/ATP
gamma2=exp(-8);
gamma=gamma0*gamma1;
gamma_t=exp(-8);
gamma_tf=exp(-8);
gamma_tb=exp(-8);
gamma_s=gamma;
K_DT=1;%relative affinity of ADP and ATP binding to KaiAKaiC
alpha(1,1:6)=1/(1+gamma1*K_DT);
KaiB_strgth=0.0;

T=[298,308];
% T = [303,313];
Efac=6;
Ep=[Efac,Efac,Efac,Efac,Efac-1,Efac+20]*8.31*T(1);Edp=Efac*8.31*T(1);Eb0=Efac*8.31*T(1);Ef6=Efac*8.31*T(1);EAf=Efac*8.31*T(1);EAb=Efac*8.31*T(1);
kp=kp_0.*exp(Ep/(8.31*T(1))).*exp(-Ep/(8.31*T(2)));
kdp=kdp_0*exp(Edp/(8.31*T(1)))*exp(-Edp/(8.31*T(2)));
b0=b0_0*exp(Eb0/(8.31*T(1)))*exp(-Eb0/(8.31*T(2)));
f6=f6_0*exp(Ef6/(8.31*T(1)))*exp(-Ef6/(8.31*T(2)));
kf=kf_0*exp(EAf/(8.31*T(1)))*exp(-EAf/(8.31*T(2)));
kb=kb_0*exp(EAb/(8.31*T(1)))*exp(-EAb/(8.31*T(2)));

incep=zeros(20,1);
incep(1)=Ct;
options=odeset('RelTol',1e-5,'NonNegative',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]);
N=19;
hwait=waitbar(0,'Waiting');
for kk=1:N
    waitbar(kk/N,hwait,'Waiting');
    added(1,kk)=-0.18+0.02*kk;
    kp6(1,kk) = added(1,kk)+kp_0(6);
    kp_vary = kp_0;
%     kp_vary(6) = kp6(1,kk);

        [t1,C1]=ode15s(@(t,C) vanZon_simple_s(C,kp_vary,kdp_0,kf_0,kb_0,f6_0,b0_0,At,Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0:0.1:500],incep,options);
        [t1,C1]=ode15s(@(t,C) vanZon_simple_s(C,kp_vary,kdp_0,kf_0,kb_0,f6_0,b0_0,At,Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0,1500],C1(end,:),options);
%         [t2,C2]=ode15s(@(t,C) vanZon_simple_s(C,kp,kdp,kf,kb,f6,b0,At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0:0.1:500],incep,options);
%         [t2,C2]=ode15s(@(t,C) vanZon_simple_s(C,kp,kdp,kf,kb,f6,b0,At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0,1500],C2(end,:),options);

    L1=floor(length(t1)/2);
%     L2=floor(length(t2)/2);

    for jj=1:length(t1)
        C1(jj,21)=Ct-sum(C1(jj,:));
    end
    
%     for jj=1:length(t2)
%         C2(jj,21)=Ct-sum(C2(jj,:));
%     end
    
    phos1=C1(:,2)+C1(:,9)+C1(:,16)+2*(C1(:,3)+C1(:,10)+C1(:,17))+3*(C1(:,4)+C1(:,11)+C1(:,18))+4*(C1(:,5)+C1(:,12)+C1(:,19))+5*(C1(:,6)+C1(:,13)+C1(:,20))+6*(C1(:,7)+C1(:,14)+C1(:,21));
    phos1=phos1/(6*Ct);
%     phos2=C2(:,2)+C2(:,9)+C2(:,16)+2*(C2(:,3)+C2(:,10)+C2(:,17))+3*(C2(:,4)+C2(:,11)+C2(:,18))+4*(C2(:,5)+C2(:,12)+C2(:,19))+5*(C2(:,6)+C2(:,13)+C2(:,20))+6*(C2(:,7)+C2(:,14)+C2(:,21));
%     phos2=phos2/(6*Ct);

%     figure(1);plot(t1(L1:end),phos1(L1:end),'-',t2(L2:end),phos2(L2:end),'-');

%     [Prd(1,kk),Upper(1,kk),Lower(1,kk)]=periods_Amp(t1(L1:end),phos1(L1:end));
%     [Prd(2,kk),Upper(2,kk),Lower(2,kk)]=periods_Amp(t2(L2:end),phos2(L2:end));
    Prd(1,kk) = periods(t1(L1:end),phos1(L1:end));
%     Prd(2,kk) = periods(t2(L1:end),phos2(L1:end));
%     if Prd(1,kk)>0
%         Q10(1,kk)=Prd(2,kk)/Prd(1,kk);
%     else
%         Prd(1,kk)=NaN;
%         Q10(1,kk)=NaN;
%     end
    
    A1=At-(C1(:,8)+C1(:,9)+C1(:,10)+C1(:,11)+C1(:,12)+C1(:,13)+C1(:,14));
%     A2=At(1,kk)-(C2(:,8)+C2(:,9)+C2(:,10)+C2(:,11)+C2(:,12)+C2(:,13)+C2(:,14));
%     %figure(2);plot(t1(L1:end),A1(L1:end));hold on;
    
    C5_total = C1(:,6) + C1(:,13);
    Csum0_4_total = C1(:,1) + C1(:,2) + C1(:,3) + C1(:,4) + C1(:,5) + C1(:,8) + C1(:,9) + C1(:,10) + C1(:,11) + C1(:,12);
    dt = 0.05;
    tt = 0:dt:t1(end-1);
    C5_total_fit = spline(t1,C5_total,tt);
    Csum0_4_total_fit = spline(t1,Csum0_4_total,tt);
    C_dephos_fit = Ct - C5_total_fit - Csum0_4_total_fit;

    [pks,ind_pks] = findpeaks(Csum0_4_total_fit);
    C5_total_thresh = C5_total_fit(ind_pks(end));
    for ii = ind_pks(end-1):ind_pks(end)
        if C5_total_fit(ii-1) > C5_total_thresh && C5_total_fit(ii) <= C5_total_thresh
            ind_trgh = ii;
            break;
        end
    end
    [pks2,ind_pks2] = findpeaks(C5_total_fit);
    ind_pks_mid = ind_pks2(ind_pks2<ind_pks(end) & ind_pks2>ind_pks(end-1));

    figure(1);
    plot(C5_total_fit(ind_trgh:ind_pks(end)),Csum0_4_total_fit(ind_trgh:ind_pks(end)),Color=[1 0 0]);
    hold on;
    plot(C5_total_fit(ind_pks(end-1):ind_pks_mid),Csum0_4_total_fit(ind_pks(end-1):ind_pks_mid),Color=[0 1 0]);
    plot(C5_total_fit(ind_pks_mid:ind_trgh),Csum0_4_total_fit(ind_pks_mid:ind_trgh),Color=[0 0 1]);
%     figure(1);
%     plot3(C5_total_fit(ind_trgh:ind_pks(end)),Csum0_4_total_fit(ind_trgh:ind_pks(end)),C_dephos_fit(ind_trgh:ind_pks(end)),Color=[1 0 0]);
%     hold on;
%     plot3(C5_total_fit(ind_pks(end-1):ind_pks_mid),Csum0_4_total_fit(ind_pks(end-1):ind_pks_mid),C_dephos_fit(ind_pks(end-1):ind_pks_mid),Color=[0 1 0]);
%     plot3(C5_total_fit(ind_pks_mid:ind_trgh),Csum0_4_total_fit(ind_pks_mid:ind_trgh),C_dephos_fit(ind_pks_mid:ind_trgh),Color=[0 0 1]);

    dX_fast1 = diff(C5_total_fit(ind_pks(end-1):ind_pks_mid));
    dY_fast1 = diff(Csum0_4_total_fit(ind_pks(end-1):ind_pks_mid));
    amp_fast1(kk) = sum(sqrt(dX_fast1.^2+dY_fast1.^2));
    tau_fast1(kk) = tt(ind_pks_mid) - tt(ind_pks(end-1));
    v_fast1(kk) = amp_fast1(kk)/tau_fast1(kk);

    dX_fast2 = diff(C5_total_fit(ind_pks_mid:ind_trgh));
    dY_fast2 = diff(Csum0_4_total_fit(ind_pks_mid:ind_trgh));
    amp_fast2(kk) = sum(sqrt(dX_fast2.^2+dY_fast2.^2));
    tau_fast2(kk) = tt(ind_trgh) - tt(ind_pks_mid);
    v_fast2(kk) = amp_fast2(kk)/tau_fast2(kk);

    amp_fast(kk) = amp_fast1(kk) + amp_fast2(kk);
    tau_fast(kk) = tt(ind_trgh) - tt(ind_pks(end-1));
    v_fast(kk) = amp_fast(kk)/tau_fast(kk);

    dX_slow = diff(C5_total_fit(ind_trgh:ind_pks(end)));
    dY_slow = diff(Csum0_4_total_fit(ind_trgh:ind_pks(end)));
    amp_slow(kk) = sum(sqrt(dX_slow.^2+dY_slow.^2));
    tau_slow(kk) = tt(ind_pks(end)) - tt(ind_trgh);
    v_slow(kk) = amp_slow(kk)/tau_slow(kk);

    dX_cycle = diff(C5_total_fit(ind_pks(end-1):ind_pks(end)));
    dY_cycle = diff(Csum0_4_total_fit(ind_pks(end-1):ind_pks(end)));
    dZ_cycle = diff(C_dephos_fit(ind_pks(end-1):ind_pks(end)));
%     v_abs = sqrt(dX_cycle.^2+dY_cycle.^2+dZ_cycle.^2)/dt;
    v_abs = sqrt(dX_cycle.^2+dY_cycle.^2)/dt;
    amp_cycle(kk) = sum(sqrt(dX_cycle.^2+dY_cycle.^2+dZ_cycle.^2));
    tau_cycle(kk) = tt(ind_pks(end)) - tt(ind_pks(end-1));
    v_cycle(kk) = amp_cycle(kk)/tau_cycle(kk);

end
toc
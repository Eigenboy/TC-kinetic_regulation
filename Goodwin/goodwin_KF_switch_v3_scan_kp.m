% Simplifed Zhou-Kim model; varying k_p
% Using central difference method and ode89 for higher accuracy of period
% sensitivities

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
k(1,1:4)= [10,10,10,10];

M = 201;

[P1,P2,Q10,DPDK,k_m,k_p] = deal(zeros(1,M));
RDP = zeros(8,M);

hwait=waitbar(0,'waiting');
for jj = 1:M
    waitbar(jj/M,hwait,'waiting');
    k_p(1,jj) = 0.885*10^(0.01*(jj-1));
    %kp onset is ~0.885;
    % k_p(1,jj) = 2.0 + 0.02*(jj-1);
    k(1,2) = k_p(1,jj);
    k(2,1:4)=k(1,1:4).*exp(Ek/(8.31*T(1))).*exp(-Ek/(8.31*T(2)));

    [~,X1]=ode89(@(t,x) Goodwin_KF_switch_v3(x,k(1,:),d(1,:),A,Kd,Km,Km2),[0,800],[1,0,0,0]);
    [t1,X1]=ode89(@(t,x) Goodwin_KF_switch_v3(x,k(1,:),d(1,:),A,Kd,Km,Km2),[0,1200],X1(end,:));
% figure(1);plot(t1,X1);
% figure(2);plot(X1(:,1),X1(:,4));hold on;
    P1(1,jj)=periods(t1(floor(length(t1)/2):end,1),X1(floor(length(t1)/2):end,4));

    if P1(1,jj)~=0
        kv=zeros(8);
        amp=0.02*[1,-1];
        for kk=1:8
            kv(1,1:4)=k(1,1:4);
            kv(1,5:8)=d(1,1:4);
            kv(1,kk)=(1+amp(1))*kv(1,kk);
            kv(2,1:4)=k(1,1:4);
            kv(2,5:8)=d(1,1:4);
            kv(2,kk)=(1+amp(2))*kv(2,kk);
       
            [~,X3]=ode89(@(t,x) Goodwin_KF_switch_v3(x,kv(1,1:4),kv(1,5:8),A,Kd,Km,Km2),[0,800],[1,0,0,0]);
            [t3,X3]=ode89(@(t,x) Goodwin_KF_switch_v3(x,kv(1,1:4),kv(1,5:8),A,Kd,Km,Km2),[0,1200],X3(end,:));
            P3=periods(t3(floor(length(t3)/2):end,1),X3(floor(length(t3)/2):end,4));
            [~,X4]=ode89(@(t,x) Goodwin_KF_switch_v3(x,kv(2,1:4),kv(2,5:8),A,Kd,Km,Km2),[0,800],[1,0,0,0]);
            [t4,X4]=ode89(@(t,x) Goodwin_KF_switch_v3(x,kv(2,1:4),kv(2,5:8),A,Kd,Km,Km2),[0,1200],X4(end,:));
            P4=periods(t4(floor(length(t4)/2):end,1),X4(floor(length(t4)/2):end,4));
            % figure(1);plot(t3,X3);
            % figure(2);plot(X3(:,1),X3(:,2));hold on;
            if P3 > 0 && P4 > 0
                RDP(kk,jj) = (P3 - P4)/(P1(1,jj)*(amp(1)-amp(2)));
            elseif P3 > 0
                RDP(kk,jj) = 1/amp(1)*(P3/P1(1,jj)-1);
            elseif P4 > 0
                RDP(kk,jj) = 1/amp(2)*(P4/P1(1,jj)-1);
            else
                RDP(kk,jj) = NaN;
            end
            clear t3 z3;
        end
        DPDK(1,jj)=sum(RDP(:,jj));
    else
        RDP(:,jj)=NaN;
        DPDK(1,jj)=NaN;
    end

    [~,X2]=ode89(@(t,x) Goodwin_KF_switch_v3(x,k(2,:),d(2,:),A,Kd,Km,Km2),[0,800],[1,0,0,0]);
    [t2,X2]=ode89(@(t,x) Goodwin_KF_switch_v3(x,k(2,:),d(2,:),A,Kd,Km,Km2),[0,1200],X2(end,:));

    P2(1,jj)=periods(t2(floor(length(t2)/2):end,1),X2(floor(length(t2)/2):end,4));
    Q10(1,jj)=P2(1,jj)/P1(1,jj);
end
toc
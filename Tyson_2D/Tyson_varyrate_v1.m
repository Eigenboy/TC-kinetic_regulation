% Tyson-2D model; varying k_1 (phosphorylation rate)

clear;
tic
v_m = 1*1.0;
k_m = 0.1*1.0;
v_p = 0.5*1.0;
k_1 = 10*0.72;
k_2 = 0.03*1.0;
k_3 = 0.1*1.0;
J = 0.05;
K = 200*1;
A = 0.1;

T=[298;308];
E0=[15,15,15,25,15,15]*8.31*T(1);
k0 = [v_m,k_m,v_p,k_1,k_2,k_3];
amp = 0.01;
P_temp = 24;

M = 200;
[k_1_vary,P1,P2,Q10,RDP,DPDK] = deal(zeros(1,M));

hwait=waitbar(0,'waiting');
for ii = 1:M
    waitbar(ii/M,hwait,'waiting');
    k = k0;
    k(2,:) = k0 .* exp(E0(1,:)/(8.31*T(1))) .* exp(-E0(1,:)/(8.31*T(2)));
    k(:,4) = k(:,4) * exp(0.01*(ii-1));
    k_1_vary(1,ii) = k(1,4);
    
    [t1,z1]=ode23s(@(t,z) Tyson_per_2D(z,k(1,1),k(1,2),k(1,3),k(1,4),k(1,5),k(1,6),J,K,A),[0,30*P_temp],[1,1]);
    P_temp1 = periods(t1,z1(:,1));
    if P_temp1 > 0
        [t1,z1]=ode23s(@(t,z) Tyson_per_2D(z,k(1,1),k(1,2),k(1,3),k(1,4),k(1,5),k(1,6),J,K,A),[0,60*P_temp1],z1(end,:));
        P_temp = periods(t1,z1(:,1));
    else
        P_temp = 0;
    end
    if P_temp > 0
        P1(1,ii) = P_temp;
    else
        break;
    end

%     figure(1);
%     plot(z1(:,1),z1(:,2));

    [~,z2]=ode23s(@(t,z) Tyson_per_2D(z,k(2,1),k(2,2),k(2,3),k(2,4),k(2,5),k(2,6),J,K,A),[0,30*P_temp],[1,1]);
    [t2,z2]=ode23s(@(t,z) Tyson_per_2D(z,k(2,1),k(2,2),k(2,3),k(2,4),k(2,5),k(2,6),J,K,A),[0,60*P_temp],z2(end,:));
    P2(1,ii) = periods(t2,z2(:,1));
    Q10(1,ii) = P2(1,ii)/P1(1,ii);
    if P1(1,ii) > 0
        for jj = 1: length(k)
            k_vary = k(1,:);
            k_vary(1,jj) = (1+amp) * k(1,jj);
            [~,z3]=ode23s(@(t,z) Tyson_per_2D(z,k_vary(1,1),k_vary(1,2),k_vary(1,3),k_vary(1,4),k_vary(1,5),k_vary(1,6),J,K,A),[0,30*P1(1,ii)],[1,1]);
            [t3,z3]=ode23s(@(t,z) Tyson_per_2D(z,k_vary(1,1),k_vary(1,2),k_vary(1,3),k_vary(1,4),k_vary(1,5),k_vary(1,6),J,K,A),[0,60*P1(1,ii)],z3(end,:));
            P3 = periods(t3,z3(:,1));
            
            % period sensitivities
            RDP(jj,ii) = 1/amp * (P3/P1(1,ii)-1);
        end
        DPDK(1,ii) = sum(RDP(:,ii));
    end
end
save('data_varyraye_v1.mat','k_1_vary','P1','P2','Q10','RDP','DPDK','v_m','k_m','v_p','k_2','k_3','J','K','A','T','E0','k0','amp');

toc

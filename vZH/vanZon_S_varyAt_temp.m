% van Zon - Hatakeyama model; varying total KaiA (At)

clear;
tic
fctr=0.333;
kp_0(1:6)=3.1*fctr;
kp_0(5)=1.4*kp_0(5);
kp_0(6)=1.2*kp_0(6)*1;
kdp_0(1:6)=9.5*fctr;
b0_0=1*fctr;
f6_0=1*fctr;
kf_0(1:7)=5*1e7*fctr;
kb_0=zeros(1,6);
strgth=10;
for jj=1:7
    kb_0(jj)=strgth^7*strgth^(jj-7)*fctr;
end
Ct=11.4;

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
Ep=[Efac,Efac,Efac,Efac,Efac,Efac+20]*8.31*T(1);Edp=Efac*8.31*T(1);Eb0=Efac*8.31*T(1);Ef6=Efac*8.31*T(1);EAf=Efac*8.31*T(1);EAb=Efac*8.31*T(1);
kp=kp_0.*exp(Ep/(8.31*T(1))).*exp(-Ep/(8.31*T(2)));
kdp=kdp_0*exp(Edp/(8.31*T(1)))*exp(-Edp/(8.31*T(2)));
b0=b0_0*exp(Eb0/(8.31*T(1)))*exp(-Eb0/(8.31*T(2)));
f6=f6_0*exp(Ef6/(8.31*T(1)))*exp(-Ef6/(8.31*T(2)));
kf=kf_0*exp(EAf/(8.31*T(1)))*exp(-EAf/(8.31*T(2)));
kb=kb_0*exp(EAb/(8.31*T(1)))*exp(-EAb/(8.31*T(2)));

incep=zeros(20,1);
incep(1)=Ct;
options=odeset('RelTol',1e-5,'NonNegative',[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]);
N=200;
hwait=waitbar(0,'Waiting');
for kk=1:N
    waitbar(kk/N,hwait,'Waiting');
    % At(1,kk)=0.6+kk*0.02;
      At(1,kk)=0.54+kk*0.02;
        [t1,C1]=ode15s(@(t,C) vanZon_simple_s(C,kp_0,kdp_0,kf_0,kb_0,f6_0,b0_0,At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0:0.1:500],incep,options);
        [t1,C1]=ode15s(@(t,C) vanZon_simple_s(C,kp_0,kdp_0,kf_0,kb_0,f6_0,b0_0,At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0,1500],C1(end,:),options);
        
        [t2,C2]=ode15s(@(t,C) vanZon_simple_s(C,kp,kdp,kf,kb,f6,b0,At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0:0.1:500],incep,options);
        [t2,C2]=ode15s(@(t,C) vanZon_simple_s(C,kp,kdp,kf,kb,f6,b0,At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0,1500],C2(end,:),options);

    L1=floor(length(t1)/2);
    L2=floor(length(t2)/2);

    for jj=1:length(t1)
        C1(jj,21)=Ct-sum(C1(jj,:));
    end
    
    for jj=1:length(t2)
        C2(jj,21)=Ct-sum(C2(jj,:));
    end
    
    phos1=C1(:,2)+C1(:,9)+C1(:,16)+2*(C1(:,3)+C1(:,10)+C1(:,17))+3*(C1(:,4)+C1(:,11)+C1(:,18))+4*(C1(:,5)+C1(:,12)+C1(:,19))+5*(C1(:,6)+C1(:,13)+C1(:,20))+6*(C1(:,7)+C1(:,14)+C1(:,21));
    phos1=phos1/(6*Ct);
    phos2=C2(:,2)+C2(:,9)+C2(:,16)+2*(C2(:,3)+C2(:,10)+C2(:,17))+3*(C2(:,4)+C2(:,11)+C2(:,18))+4*(C2(:,5)+C2(:,12)+C2(:,19))+5*(C2(:,6)+C2(:,13)+C2(:,20))+6*(C2(:,7)+C2(:,14)+C2(:,21));
    phos2=phos2/(6*Ct);

    figure(1);plot(t1(L1:end),phos1(L1:end),'-',t2(L2:end),phos2(L2:end),'-');

    [Prd(1,kk),Upper(1,kk),Lower(1,kk)]=periods_Amp(t1(L1:end),phos1(L1:end));
    [Prd(2,kk),Upper(2,kk),Lower(2,kk)]=periods_Amp(t2(L2:end),phos2(L2:end));
    if Prd(1,kk)>0
        Q10(1,kk)=Prd(2,kk)/Prd(1,kk);
    else
        Prd(1,kk)=NaN;
        Q10(1,kk)=NaN;
    end
    
    A1=At(1,kk)-(C1(:,8)+C1(:,9)+C1(:,10)+C1(:,11)+C1(:,12)+C1(:,13)+C1(:,14));
    A2=At(1,kk)-(C2(:,8)+C2(:,9)+C2(:,10)+C2(:,11)+C2(:,12)+C2(:,13)+C2(:,14));
    %figure(2);plot(t1(L1:end),A1(L1:end));hold on;

    %varyrates:
    if Prd(1,kk) > 0
        amp=0.05;
        kv=zeros(28);
        for ii=1:28
            kv(ii,1:6)=kp_0;
            kv(ii,7:12)=kdp_0;
            kv(ii,13:19)=kf_0;
            kv(ii,20:26)=kb_0;
            kv(ii,27)=f6_0;
            kv(ii,28)=b0_0;
            kv(ii,ii)=(1+amp)*kv(ii,ii);
            [t3,C3]=ode15s(@(t,C) vanZon_simple_s(C,kv(ii,1:6),kv(ii,7:12),kv(ii,13:19),kv(ii,20:26),kv(ii,27),kv(ii,28),At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s),[0:0.1:500],incep,options);
            [t3,C3]=ode15s(@(t,C) vanZon_simple_s(C,kv(ii,1:6),kv(ii,7:12),kv(ii,13:19),kv(ii,20:26),kv(ii,27),kv(ii,28),At(1,kk),Ct,gamma,gamma2,gamma_tf,gamma_tb,alpha,gamma_s/(1+amp)),[0,1500],C3(end,:),options);
            for jj=1:length(t3)
                C3(jj,21)=Ct-sum(C3(jj,:));
            end
        
            phos3=C3(:,2)+C3(:,9)+C3(:,16)+2*(C3(:,3)+C3(:,10)+C3(:,17))+3*(C3(:,4)+C3(:,11)+C3(:,18))+4*(C3(:,5)+C3(:,12)+C3(:,19))+5*(C3(:,6)+C3(:,13)+C3(:,20))+6*(C3(:,7)+C3(:,14)+C3(:,21));
            phos3=phos3/(6*Ct);
            Prd3=periods(t3,phos3);

            % period sensitivities
            RDP(ii,kk)=(Prd3/Prd(1,kk)-1)/amp;
            clear t3 C3;
        end
        DPDK(1,kk) = sum(RDP(:,kk));
    else
        RDP(1:28,kk) = NaN;
        DPDK(1,kk) = NaN;
    end
    
end
toc

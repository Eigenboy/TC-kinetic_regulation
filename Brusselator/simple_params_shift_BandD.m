% Brusselator; varying both [B] and [D], and computing Q10, energy dissipation, and
% period sensitivity of k_2

clear;
tic
a=1;
Mt=1;
T=[298;308];

k0=[1,2,1,1,10^(-6),2*10^(-6),10^(-6),10^(-6)];
E0=[15,25,20,20,20,20,20,20]*8.31*T(1);

N = 73;
M = 73;
[B,D,P1,Dsp1,P2,Q10,K,RDP2,X_min,X_max,Y_min,Y_max]=deal(zeros(N,M));
hwait=waitbar(0,'waiting');

for jj = 1:N
    waitbar(jj/N,hwait,'waiting');
    for kk=1:M
        B(jj,kk)=exp(- 0.08*jj + 0.1);
        D(jj,kk)=exp(0.25*kk - 3.25);
    
        k=zeros(2,8);
        for ii=1:8
            if ii==2
               k(1,ii)=k0(1,ii)*B(jj,kk);
               k(2,ii)=k0(1,ii)*exp(E0(1,ii)/(8.31*T(1)))*exp(-E0(1,ii)/(8.31*T(2)))*B(jj,kk);
            elseif ii==6
               k(1,ii)=k0(1,ii)*D(jj,kk);
               k(2,ii)=k0(1,ii)*exp(E0(1,ii)/(8.31*T(1)))*exp(-E0(1,ii)/(8.31*T(2)))*D(jj,kk);
            else
               k(1,ii)=k0(1,ii);
               k(2,ii)=k0(1,ii)*exp(E0(1,ii)/(8.31*T(1)))*exp(-E0(1,ii)/(8.31*T(2)));
            end
        end
    
        x_st=(k(:,1)+k(:,8))*a./(k(:,5)+k(:,4));
        y_st=(k(:,2)+k(:,7).*x_st(:,1).^2).*x_st(:,1)./(k(:,6)+k(:,3).*x_st(:,1).^2);
        lambda1=YacobiEigen1(k(1,:),x_st(1,1),y_st(1,1));
        lambda2=YacobiEigen1(k(2,:),x_st(2,1),y_st(2,1));
        if lambda1(1,1)>0
%         if abs(imag(lambda1(1,1)))>1e-10
%             omega1=abs(imag(lambda1(1,1))); 
%         else
%             omega1=8*pi*(k(1,1)+k(1,8))^2*k(1,3)/k(1,2)^2;
%         end
           if kk == 1
               ksai=k(1,6)*(k(1,2) + k(1,4) + k(1,5)).^2/(k(1,3)*(k(1,1)+k(1,8))^2);
               omega1=4*pi * k(1,6) / (sqrt(1+ksai)-1);
           elseif P1(jj,kk-1) > 0
               omega1 = 2*pi/P1(jj,kk-1);
           else
               ksai=k(1,6)*(k(1,2) + k(1,4) + k(1,5)).^2/(k(1,3)*(k(1,1)+k(1,8))^2);
               omega1=4*pi * k(1,6) / (sqrt(1+ksai)-1);
           end
           [t1,z1]=ode23s(@(t,x) bruss_new(x,k(1,:),a),[0,40*pi/omega1],[1.1*x_st(1,1),1.1*y_st(1,1)]);
           [t1,z1]=ode23s(@(t,x) bruss_new(x,k(1,:),a),[0,80*pi/omega1],z1(end,:));

           [P1(jj,kk),Dsp1(jj,kk),Diss1,i1]=prd_and_dsp_2reactions(t1,z1,k(1,:));
           X_min(jj,kk)=min(z1(:,1));Y_min(jj,kk)=min(z1(:,2));
           X_max(jj,kk)=max(z1(:,1));Y_max(jj,kk)=max(z1(:,2));
           
            % period sensitivity of k_2
           if P1(jj,kk)~=0
                amp=0.01;
                kv2=k(1,:);
                kv2(1,2)=(1+amp)*kv2(1,2);
                [t3,z3]=ode23s(@(t,x) bruss_new(x,kv2,a),[0,40*pi/omega1],[1.1*x_st(1,1),1.1*y_st(1,1)]);
                [t3,z3]=ode23s(@(t,x) bruss_new(x,kv2,a),[0,80*pi/omega1],z3(end,:));
                P3=periods(t3,z3(:,1));
                RDP2(jj,kk)=1/amp*(P3/P1(jj,kk)-1);
                clear t3 z3;
           else
                RDP2(jj,kk)=NaN;
           end
           
           [t2,z2]=ode23s(@(t,x) bruss_new(x,k(2,:),a),[0,40*pi/omega1],[1.1*x_st(1,1),1.1*y_st(1,1)]);
           [t2,z2]=ode23s(@(t,x) bruss_new(x,k(2,:),a),[0,80*pi/omega1],z2(end,:));
           [P2(jj,kk),Dsp2,Diss2,i2]=prd_and_dsp_2reactions(t2,z2,k(2,:));
           Q10(jj,kk)=P2(jj,kk)/P1(jj,kk);
            
        else
           P1(jj,kk)=NaN;
           Dsp1(jj,kk)=NaN;
           P2(jj,kk)=NaN;
           Q10(jj,kk)=NaN;
           RDP2(jj,kk)=NaN;
        end
    end
end
toc
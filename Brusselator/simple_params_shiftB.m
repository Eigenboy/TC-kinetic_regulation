% Brusselator; varying [B] and computing Q10 and
% period sensitivities
% all reverse reactions at set to be extremely small

clear;
tic
a=1;
D=1;
T=[298;308];

% rates at 298K
k0=[1,2,1,1,10^(-8),2*10^(-8),10^(-8),10^(-8)];
% activation energies
E0=[15,25,20,20,20,20,20,20]*8.31*T(1);

M=93;
[P1,Dsp1,P2,Dsp2,Q10,DPDK,Z_min,Z_max]=deal(zeros(1,M));
RDP=zeros(8,M);
hwait=waitbar(0,'waiting');

 for kk=1:M
    waitbar(kk/M,hwait,'waiting');
    B(1,kk)=exp(0.05*kk-0.03);

    k=zeros(2,8);
    for ii=1:8
        if ii==2
            k(1,ii)=k0(1,ii)*B(1,kk);
            k(2,ii)=k0(1,ii)*exp(E0(1,ii)/(8.31*T(1)))*exp(-E0(1,ii)/(8.31*T(2)))*B(1,kk);
        elseif ii==6
           k(1,ii)=k0(1,ii)*D;
           k(2,ii)=k0(1,ii)*exp(E0(1,ii)/(8.31*T(1)))*exp(-E0(1,ii)/(8.31*T(2)))*D;
        else
           k(1,ii)=k0(1,ii);
           k(2,ii)=k0(1,ii)*exp(E0(1,ii)/(8.31*T(1)))*exp(-E0(1,ii)/(8.31*T(2)));
        end
    end
    
    % eigenvalues at the fixed point
    x_st=(k(:,1)+k(:,8))*a./(k(:,5)+k(:,4));
    y_st=(k(:,2)+k(:,7).*x_st(:,1).^2).*x_st(:,1)./(k(:,6)+k(:,3).*x_st(:,1).^2);
    lambda1=YacobiEigen1(k(1,:),x_st(1,1),y_st(1,1));
    lambda2=YacobiEigen1(k(2,:),x_st(2,1),y_st(2,1));
    
    if lambda1(1,1)>0
        % estimation of the angular velocity
        if abs(imag(lambda1(1,1)))>1e-10
            omega1=abs(imag(lambda1(1,1))); 
        else
            omega1=8*pi*(k(1,1)+k(1,8))^2*k(1,3)/k(1,2)^2;
        end
       [t1,z1]=ode23s(@(t,x) bruss_new(x,k(1,:),a),[0,30*pi/omega1],[1.1*x_st(1,1),1.1*y_st(1,1)]);
       [t1,z1]=ode23s(@(t,x) bruss_new(x,k(1,:),a),[0,100*pi/omega1],z1(end,:));
       %plot(z1(:,1),z1(:,2));
       [P1(1,kk),Dsp1(1,kk),diss1,i1]=prd_and_dsp_S(t1,z1,k(1,:),a);
       Z_min(1,kk)=min(z1(:,1));Z_min(2,kk)=min(z1(:,2));
       Z_max(1,kk)=max(z1(:,1));Z_max(2,kk)=max(z1(:,2));
       
       % computing period sensitivities
       if P1(1,kk)~=0
            kv=zeros(8);
            amp=0.01;
            for ii=1:8
                kv(ii,:)=k(1,:);
                kv(ii,ii)=(1+amp)*kv(ii,ii);
            end

            for ii=1:8
               [t3,z3]=ode23s(@(t,x) bruss_new(x,kv(ii,:),a),[0,30*pi/omega1],[1.1*x_st(1,1),1.1*y_st(1,1)]);
               [t3,z3]=ode23s(@(t,x) bruss_new(x,kv(ii,:),a),[0,100*pi/omega1],z3(end,:));
               P3=periods(t3,z3(:,1));

               % period sensitivities
               RDP(ii,kk)=1/amp*(P3/P1(1,kk)-1);
               clear t3 z3;
            end
            % verifying the sum rule
            DPDK(1,kk)=sum(RDP(:,kk));
        else
            RDP(:,kk)=NaN;
            DPDK(1,kk)=NaN;
       end
       
       [t2,z2]=ode23s(@(t,x) bruss_new(x,k(2,:),a),[0,30*pi/omega1],[1.1*x_st(1,1),1.1*y_st(1,1)]);
       [t2,z2]=ode23s(@(t,x) bruss_new(x,k(2,:),a),[0,100*pi/omega1],z2(end,:));
       [P2(1,kk),Dsp2(1,kk),diss2,i2]=prd_and_dsp_S(t2,z2,k(2,:),a);
       Q10(1,kk)=P2(1,kk)/P1(1,kk);
       clear t2 z2;
        
    else
       P1(1,kk)=NaN;
       Dsp1(1,kk)=NaN;
       P2(1,kk)=NaN;
       Dsp2(1,kk)=NaN;
       Q10(1,kk)=NaN;
       RDP(:,kk)=NaN;
       DPDK(1,kk)=NaN;
    end
end
toc
function [avT,W,Diss,i1]=prd_and_dsp_2reactions(t,z,k)
% compute both period and the energy dissipation for brusselator
% only dissipation of the reactions between X and Y is considered

N=length(t);
jj=1;
T(jj)=t(1);
ma=max(z(:,1));
mi=min(z(:,1));
mid=0.5*(ma+mi);  
i1=0;

%count the events of passing the middle of z(:,1).
for ii=1:N-1
    if (z(ii,1)<mid)&&(z(ii+1,1)>mid)
        T(jj)=t(ii);
        i1(jj)=ii;
        jj=jj+1;
    end
end
M=length(T);

diss = zeros(i1(M)-i1(1)+1,2);
if (M>5)&&((ma-mi)/mean(z(:,1))>0.03)
    for ii=i1(1):i1(M)
        Jf=zeros(2,1);                   Jb=zeros(2,1);
        Jf(1)=k(2)*z(ii,1);              Jb(1)=k(6)*z(ii,2);
        Jf(2)=k(3)*z(ii,1)^2*z(ii,2);    Jb(2)=k(7)*z(ii,1)^3;
        diss(ii-i1(1)+1,1)= tdiss(Jf(1),Jb(1));
        diss(ii-i1(1)+1,2)= tdiss(Jf(2),Jb(2));
    end
    Diss(1,1)=trapz(t(i1(1):i1(M)),diss(:,1))/(M-1);
    Diss(2,1)=trapz(t(i1(1):i1(M)),diss(:,2))/(M-1);
    avT=(T(M)-T(1))/(M-1);  
    W=sum(Diss);    
else
    avT=0;
    Diss=[0;0];
    W=0;
end
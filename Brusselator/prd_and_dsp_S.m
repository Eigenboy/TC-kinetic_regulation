function [avT,W,diss,i1]=prd_and_dsp_S(t,z,k,a)
%period and dissipation for brusselator
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

if (M>5)&&((ma-mi)/mean(z(:,1))>0.03)
    for ii=i1(1):i1(M)
        Jf=zeros(4,1);                   Jb=zeros(4,1);
        Jf(1)=k(1)*a;                    Jb(1)=k(5)*z(ii,1);
        Jf(2)=k(2)*z(ii,1);              Jb(2)=k(6)*z(ii,2);
        Jf(3)=k(3)*z(ii,1)^2*z(ii,2);    Jb(3)=k(7)*z(ii,1)^3;
        Jf(4)=k(4)*z(ii,1);              Jb(4)=k(8)*a;
        d(ii-i1(1)+1)=sum((Jf(Jf>0)-Jb(Jf>0)).*log(Jf(Jf>0)./Jb(Jf>0)));
    end
    diss=d';
    Diss=trapz(t(i1(1):i1(M)),d);
    avT=(T(M)-T(1))/(M-1);  
    W=Diss/(M-1);    
else
    avT=0;
    W=0;
end
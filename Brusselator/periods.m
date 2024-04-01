function avT=periods(t,z)
%period of one oscillatory array
N=length(t);
jj=1;
T(jj)=t(1);
mid=0.5*(max(z)+min(z));
for ii=1:N-1
    if (z(ii)<mid)&&(z(ii+1)>mid)
        T(jj)=t(ii);
        jj=jj+1;
    end
end
M=length(T);

if (M>5)&&((max(z(:,1))-min(z(:,1)))/mean(z(:,1))>0.03)
    avT=(T(M)-T(1))/(M-1);
else
    avT=0;
end
        
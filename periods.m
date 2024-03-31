function avT=periods(t,z)
% computing the period of a periodic sequence
% t is the time sequence and z is the corresponding periodic sequence
% recording the events when z ascends across the middle value

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
if abs((max(z)-min(z))/mean(z))>0.10 && M>3
    avT=(T(M)-T(1))/(M-1);
else
    avT=0;
end
        
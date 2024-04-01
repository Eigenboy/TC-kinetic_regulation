function avT=periods(t,phos)
% computing the period of a periodic sequence
% t is the time sequence and z is the corresponding periodic sequence
% recording the events when z ascends across the middle value

N=length(t);
jj=1;
T(jj)=t(1);
mid=0.5*(max(phos)+min(phos));
for ii=1:N-1
    if (phos(ii)<mid)&&(phos(ii+1)>mid)
        T(jj)=t(ii);
        jj=jj+1;
    end
end
M=length(T);
if (M>4)&&((max(phos(:,1))-min(phos(:,1)))/mean(phos(:,1))>0.002)
    avT=(T(M)-T(1))/(M-1);
else
    avT=0;
end
        
function [avT,mx,mn]=periods_Amp(t,phos)
% computing the period and the amplitude of a periodic sequence
% t is the time sequence and z is the corresponding periodic sequence
% recording the events when z ascends across the middle value
% recording peaks and troughs by findpeaks

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
[peaks,t_peaks]=findpeaks(phos,t);
[troughs,t_troughs]=findpeaks(-phos,t);
M=length(T);
if (M>4)&&((max(phos(:,1))-min(phos(:,1)))/mean(phos(:,1))>0.002)
    avT=(T(M)-T(1))/(M-1);
    troughs=-troughs;
    mx=mean(peaks);
    mn=mean(troughs);
else
    avT=0;
    mx=mean(phos);
    mn=mx;
end

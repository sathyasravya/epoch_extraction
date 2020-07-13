clc;
clear all;
close all;
[x, fs] = audioread("09a07Wb.wav"); %reads data file
[egg, fs1] = audioread("09a07Wbxx.wav");%reads egg wave file
N = length(x);% length of the signal
signal_duration = N/fs;% in seconds
x_shifted = [0;x];%shifted to right by 1
diffSig_x = x-x_shifted(1:N,1);% x(n)-x(n-1)
del_F=10;% frequency spacing, del F
K = fs/(2*del_F);% No of frequency components
y = zeros(K,N);% Matrix to store filtered signal
r=0.995;% r value for filter
%% performs filtering
for k=1:K
    fk = k*del_F;
    fk_bar = fs/2 - fk;
   for n=2:N
        mfactor = exp(-1i*2*pi*fk_bar*(n)/fs);% vector
        xkn = diffSig_x(n)*mfactor;
        y(k,n) = -r*y(k,n-1)+xkn;
    end
end
v=abs(y);
%% Spectral mean and spectral variance
sumvec = sum(v,1);
meanvec = sumvec'/K;
normvec = zeros(K,N);
for it = 1:N
    for k1 = 1:K
        normvec(k1,it) = v(k1,it)/sumvec(it);
    end
end
stdvec = std(normvec,0,1)';
%% Epoch Extraction using ZFF, mZFF(modified ZFF)
[Wlen]=xcorrWinLen(x,fs); % generates a window length
[zfSig]=zeroFreqFilter(x,fs,Wlen); % ZFF
[gci]=epochExtract(zfSig);% extracts epochs from ZFF signal
%% Intervals of epoch presence derived from mZFF signal(zero crossings)
Window=zeros(size(x));% Rectangular window from zero crossings of mZFF
No_epochs=find(gci==1);% No of epoch locations
for i=1:length(No_epochs)
    index=No_epochs(i);
    Window(index-8:index+8)=1;
end
%% Generates epoch locations for spectral mean and spectral variance of SFF 
mZFF_SFF_mean=zeros(size(x));
mZFF_SFF_variance=zeros(size(x));
for j=1:length(No_epochs)
    index=No_epochs(j);
    mean_pdct1=Window(index-8:index+8).*meanvec(index-8:index+8);% product of window and spectral mean
    [Max1, maxindx1]=max(mean_pdct1);% finds peak in the interval
    mZFF_SFF_mean(index-8+maxindx1-1)=Max1;
    var_pdct2=Window(index-8:index+8).*stdvec(index-8:index+8);% product of window and spectral variance
    [Max2, maxindx2]=max(var_pdct2);
    mZFF_SFF_variance(index-8+maxindx2-1)=Max2;
end
%% Ground truth
egg=egg(:,2); % EGG signal
degg = diff(egg); degg(end+1)=degg(end);    %dEGG signal
time=(1:length(x))*1000/fs;    % in milliseconds
%% plotting
figure(1)
subplot(5,1,1)
plot(x(6001:6400,1));
xlim([0,400])
xlabel('time(ms)'); ylabel('Amplitude')
title('A Segment of voiced speech signal')
subplot(5,1,2)
plot(diffSig_x(6001:6400,1));
xlim([0,400])
xlabel('time(ms)'); ylabel('Amplitude')
title('Segment of Pre-emphasis signal')
subplot(5,1,3)
imagesc(1:400,1:16000,v(:,6001:6400));
ylabel('Frequency')
xlabel('time(ms)');
title('SFF Spectrogram')
subplot(5,1,4)
plot(stdvec(6001:6400,1))
xlim([0,400])
xlabel('time(ms)');
title('Spectral Variance')
subplot(5,1,5)
plot(degg(6001:6400,1)./max(degg))
xlim([0,400]);
xlabel('time(ms)'); ylabel('Amplitude');
title('Ground truth')

figure(2)
subplot(3,1,1)
plot(meanvec(6001:6400,1))
xlim([0,400]);
xlabel('time(ms)');
title('Spectral Gain')


subplot(3,1,2)
plot(stdvec(6001:6400,1))
xlim([0,400])
xlabel('time(ms)');
title('Spectral Variance')
subplot(3,1,3)
plot(degg(6001:6400,1)./max(degg))
xlim([0,400]);
xlabel('time(ms)'); ylabel('Amplitude');
title('Ground truth')

figure(3)
subplot(3,3,1)
plot(x(6001:6400,1)/max(x));xlim([0,400]);xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude');title('A Segment of voiced speech signal')
subplot(3,3,2)
plot(zfSig(6001:6400,1)./max(zfSig));xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude'); title('mZFF Signal');
subplot(3,3,3)
stem(gci(6001:6400),'Marker','none'); xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude'); title('Epochs from mZFF Signal'); ylim([0 1.5]);
subplot(3,3,4)
plot( Window(6001:6400));xlim([0,400]);xlim([0,400]);ylim([0,1.2]);xlabel('time(ms)'); ylabel('Amplitude'); title('Intervals of epoch presence deried from mZFF signal');
subplot(3,3,5)
plot( meanvec(6001:6400)/max(meanvec));xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude'); title('Spectral Gain, SFF(mean)');
subplot(3,3,6)
plot( mZFF_SFF_mean(6001:6400)/max(mZFF_SFF_mean));xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude'); title('mZFF+SFF(mean)');
subplot(3,3,7)
plot( stdvec(6001:6400,1)/max(stdvec));xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude'); title('Spectral Variance,SFF(variance)');
subplot(3,3,8)
plot( mZFF_SFF_variance(6001:6400)/max(mZFF_SFF_variance));xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude'); title('mZFF+SFF(var)');
subplot(3,3,9)
plot(degg(6001:6400,1)./max(degg));xlim([0,400]);xlabel('time(ms)'); ylabel('Amplitude');title('Ground truth')

figure(4)
plot(x(6001:6400,1));
xlim([0,400])
xlabel('time(ms)'); ylabel('Amplitude')
title('A Segment of voiced speech signal')
%% ====================================================== Part 1: Plotting 
clear;clc;

data    = load('data/ECG_Sig.mat');
sig     = data.Sig;
fs      = data.sfreq;
time    = (0:size(sig,1)-1)/fs;

f11 = figure;
f11.WindowState = 'maximized';

subplot(2,1,1)
plot(time,sig(:,1))
title('Signal on channel 1')
xlabel('Time (s)'); ylabel('Amplitude')
grid on
subplot(2,1,2)
plot(time,sig(:,2))
title('Signal on channel 2')
xlabel('Time (s)'); ylabel('Amplitude')
grid on

time1 = find(time<371.2 &  time> 370.2);
beat1 = sig(time1,1);
time1 = time(time1);

time2 = find(time<475.5 &  time> 474.5);
beat2 = sig(time2,1);
time2 = time(time2);

f12 = figure;
f12.WindowState = 'maximized';
plot(time1,beat1);
x = [370.5,370.59,370.635,370.75,370.87];
y = [-0.185,-0.185,0.775,-0.25,-0.055];
title('Sample heartbeat 1')
hold on; grid on
plot(x,y,'O');

f13 = figure;
f13.WindowState = 'maximized';
plot(time2,beat2);
x = [474.94,475.05,475.105,475.135,475.37];
y = [-0.135,-0.23,0.93,-0.285,0.03];
title('Sample heartbeat 2')
hold on; grid on
plot(x,y,'O');


%% ====================================================== Part 2: Plotting 
annots = data.ANNOTD;
rpeaks = data.ATRTIMED;

f2 = figure;
f2.WindowState = 'maximized';

subplot(2,1,1)
plot(time,sig(:,1))
title('Signal on channel 1')
xlabel('Time (s)'); ylabel('Amplitude')
grid on; hold on
for i = 1:length(annots)
    text(rpeaks(i),1,num2str(annots(i)))
end

subplot(2,1,2)
plot(time,sig(:,2))
title('Signal on channel 2')
xlabel('Time (s)'); ylabel('Amplitude')
grid on; hold on
for i = 1:length(annots)
    text(rpeaks(i),0,num2str(annots(i)))
end


%% ========================================================= Part 3: Types 
labels  = unique(data.ANNOTD);
figs    = cell(1,length(labels));

for i = 1 : length(labels)
    idx  = find(annots==labels(i));
    idx  = idx(1);
    rp   = rpeaks(idx);
    
    T = find(time>rp-0.3 & time<rp+0.6);
    B = sig(T,1);
    T = time(T);
    
    figs{i} = figure;
    figs{i}.WindowState = 'maximized';
    
    plot(T,B)
    title(['Sample heartbeat for label ',num2str(i)])
    grid on
    xlabel('Time (s)'); ylabel('Amplitude')
    xline(rp,'-','R');
end


%% ====================================================== Part 4: DFT,STFT 
% abnormality candiates: 555,1242,1272,1288,1301,1362
normal = 5;
abnorm = 1288;

tn = find(time>rpeaks(normal)-0.3 & time<rpeaks(normal+2)+0.35);
ta = find(time>rpeaks(abnorm)-0.15 & time<rpeaks(abnorm+2)+0.4);

sig_n = sig(tn,1);
sig_a = sig(ta,1);

tn = time(tn);
ta = time(ta);

f41 = figure;
f41.WindowState = 'maximized';

subplot(2,1,1)
plot(tn,sig_n)
grid on
xlabel('Time (s)'); ylabel('Amplitude')
title('Sample normal partition on channel 1')

subplot(2,1,2)
plot(ta,sig_a)
grid on
xlabel('Time (s)'); ylabel('Amplitude')
title('Sample abnormal partition on channel 1')

N_n = length(sig_n);
Y1 = fft(sig_a);
Z2 = abs(Y1/N_n);
Z1 = Z2(1:floor(N_n/2)+1);
Z1(2:end-1) = 2*Z1(2:end-1);
Z_n = Z1;

N_a = length(sig_a);
Y2 = fft(sig_a);
Z2 = abs(Y2/N_a);
Z1 = Z2(1:floor(N_a/2)+1);
Z1(2:end-1) = 2*Z1(2:end-1);
Z_a = Z1;


f42 = figure;
f42.WindowState = 'maximized';

subplot(2,1,1)
fr = fs*(0:floor(N_n/2))/N_n; % frequency range
plot(fr,Z_n) 
grid on
title('Single-Sided Amplitude Spectrum of Normal Partition on Channel 1')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 50])

subplot(2,1,2)
fr = fs*(0:floor(N_a/2))/N_a; % frequency range
plot(fr,Z_a) 
grid on
title('Single-Sided Amplitude Spectrum of Abnormal Partition on Channel 1')
xlabel('Frequency (Hz)')
ylabel('Amplitude')
xlim([0 50])

f43 = figure;
f43.WindowState = 'maximized';

window = floor(fs/2);
noverlap = floor(fs/4);
nfft = floor(fs/2);

subplot(1,2,1)
spectrogram(sig_n,hamming(window),noverlap,nfft,fs,'onesided');
xlim([0 50])
title('Single-Sided Spectogram of Normal Partition on Channel 1')

subplot(1,2,2)
spectrogram(sig_a,hamming(window),noverlap,nfft,fs,'onesided');
xlim([0 50])
title('Single-Sided Spectogram of abormal Partition on Channel 1')

